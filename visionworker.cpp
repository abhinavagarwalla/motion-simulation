#include "visionworker.h"
#include "robocup_ssl_client.h"
#include "timer.h"
#include <QDebug>
#include "messages_robocup_ssl_detection.pb.h"
#include "messages_robocup_ssl_geometry.pb.h"
#include "messages_robocup_ssl_wrapper.pb.h"
#include "pose.h"
#include "vision-velocity.hpp"
#include <assert.h>
#include <QtCore/QFile>
#include <iostream>
#include <stdio.h>
#include "opencv2/video/tracking.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include <QFile>
#include <QTextStream>

using namespace cv;
using namespace std;

void printRobotInfo(const SSL_DetectionRobot & robot) {
    printf("CONF=%4.2f ", robot.confidence());
    if (robot.has_robot_id()) {
        printf("ID=%3d ",robot.robot_id());
    } else {
        printf("ID=N/A ");
    }
    printf(" HEIGHT=%6.2f POS=<%9.2f,%9.2f> ",robot.height(),robot.x(),robot.y());
    if (robot.has_orientation()) {
        printf("ANGLE=%6.3f ",robot.orientation());
    } else {
        printf("ANGLE=N/A    ");
    }
    printf("RAW=<%8.2f,%8.2f>\n",robot.pixel_x(),robot.pixel_y());
}
inline void linearTransform(double &x, double &y, double &newangle)
{
  double tempx, tempy;
  tempx = (y*(-HALF_FIELD_MAXX))/(2050.0);
  tempy = ((x + 3025.0/2)*(HALF_FIELD_MAXY))/(3025.0/2.0);
  x = tempx;
  y = tempy;
  newangle = normalizeAngle(newangle+PI/2);
  //y=tempy-800;
}
VisionWorker::VisionWorker(QObject *parent) :
    QObject(parent)
{
    detectionCount = 0;
    isTeamYellow = false;
}

void VisionWorker::setup(QThread *cThread, BeliefState *bs_, QMutex *bsMutex_, bool isTeamYellow_)
{
    myThread = cThread;
    bs = bs_;
    bsMutex = bsMutex_;
    isTeamYellow = isTeamYellow_;
    connect(cThread, SIGNAL(started()), this, SLOT(onEntry()));
    connect(cThread, SIGNAL(destroyed()), this, SLOT(onExit()));
    for (int i = 0; i < MAX_BS_Q; i++) {
        bsQ.push(make_pair(BeliefState(), 0));
    }
}

void VisionWorker::onEntry()
{
    RoboCupSSLClient client;
    client.open(true);
    SSL_WrapperPacket packet;

    qDebug() << "here" << endl;
    while(true) {
        if (client.receive(packet)) {
            if (packet.has_detection()) {
                double nowTime = packet.detection().t_capture();
                double timeMs = (nowTime - bsQ.front().second)*1000.0;
                if (timeMs <= 0) {
                    timeMs = 0.001; // this can happen at the beginning. who cares.
                    qDebug() << "timeMs is negative!";
                }
                detectionCount++;
                if(detectionCount == maxDetectionCount) {
                    detectionCount = 0;
                    bsMutex->lock();
                    bs->ballIsPresent = false;
                    for(int i = 0; i < 5; i++) {
                        bs->homeIsPresent[i] = false;
                        bs->awayIsPresent[i] = false;
                    }
                    bsMutex->unlock();
                }
                SSL_DetectionFrame detection = packet.detection();
                //Display the contents of the robot detection results:
//                double t_now = GetTimeSec();

                int balls_n = detection.balls_size();
                int robots_blue_n =  detection.robots_blue_size();
                int robots_yellow_n =  detection.robots_yellow_size();
                float rhoKFx, DrhoKFx, thetaKFx, xCoordinateOfObject,newx,lastPosx = 0;
                //Ball info:
                if(balls_n) {
                    // only take 1st ball into account;
                    SSL_DetectionBall ball = detection.balls(0);
                    double ballX, ballY, dummy = 0;
                    float ballVx, ballVy;
                    ballX = ball.x();
                    ballY = ball.y();
                    linearTransform(ballX, ballY, dummy);
                    VisionVelocity::calcBallVelocity(ballX - bsQ.front().first.ballX, ballY - bsQ.front().first.ballY, timeMs, ballVx, ballVy);
                    bsMutex->lock();
                    bs->ballX = ballX;
                    bs->ballY = ballY;
                    bs->ballIsPresent = true;
                    bs->ballVx = ballVx;
                    bs->ballVy = ballVy;
                    bsMutex->unlock();
                    QString outputFilename = "/home/robocup/qt_test.txt";
                    QFile outputFile(outputFilename);
                    QTextStream outStream(&outputFile);
                    outputFile.open(QIODevice::WriteOnly | QIODevice::Append);
                    outStream << ballX <<" "<< bsQ.front().first.ballX <<" "<< ballVx <<" "<< timeMs*0.001<<"\n" ;

                    float avgvx=0,avgvy=0;

                        velxq.push_back(ballVx);
                        velyq.push_back(ballVy);
                        for(int i=1;i<4;i++)
                        {
                            if(i>velxq.size()) break;
                            avgvx+=velxq[i-1];
                            avgvy+=velyq[i-1];

                        }
                        bs->ballVx = avgvx/velxq.size();
                        bs->ballVy = avgvy/velyq.size();
                        while(velxq.size() > 3)
                            velxq.pop_front();
                        while(velyq.size() > 3)
                            velyq.pop_front();

                }
                //Blue robot info:
                // currently only handling blue robot info as home
                for (int i = 0; i < robots_blue_n; i++) {
                    double botX, botY, botTheta;
                    float botVl, botVr;
                    SSL_DetectionRobot robot = detection.robots_blue(i);
                    int botId = robot.robot_id();
                    botX = robot.x();
                    botY = robot.y();
                    if(robot.has_orientation())
                        botTheta = robot.orientation();
                    else
                        botTheta = 0;
                    linearTransform(botX, botY, botTheta);
                    if (isTeamYellow) {
                        VisionVelocity::BotPose p1(bsQ.front().first.awayX[botId], bsQ.front().first.awayY[botId], bsQ.front().first.awayTheta[botId]);
                        VisionVelocity::BotPose p2(botX, botY, botTheta);
                        VisionVelocity::calcBotVelocity(p1, p2, timeMs, botVl, botVr);
                        bsMutex->lock();
                        bs->awayX[botId] = botX;
                        bs->awayY[botId] = botY;
                        bs->awayTheta[botId] = botTheta;
                        bs->awayVl[botId] = botVl;
                        bs->awayVr[botId] = botVr;
                        bs->awayIsPresent[botId] = true;
                        bsMutex->unlock();
                    } else {
                        VisionVelocity::BotPose p1(bsQ.front().first.homeX[botId], bsQ.front().first.homeY[botId], bsQ.front().first.homeTheta[botId]);
                        VisionVelocity::BotPose p2(botX, botY, botTheta);
                        VisionVelocity::calcBotVelocity(p1, p2, timeMs, botVl, botVr);
                        bsMutex->lock();
                        bs->homeX[botId] = botX;
                        bs->homeY[botId] = botY;
                        bs->homeTheta[botId] = botTheta;
                        bs->homeVl[botId] = botVl;
                        bs->homeVr[botId] = botVr;
                        bs->homeIsPresent[botId] = true;
                        bsMutex->unlock();
                    }
                }

                //Yellow robot info:
                // currently yellow is away
                for (int i = 0; i < robots_yellow_n; i++) {
                    double botX, botY, botTheta;
                    float botVl, botVr;
                    SSL_DetectionRobot robot = detection.robots_yellow(i);
                    int botId = robot.robot_id();
                    botX = robot.x();
                    botY = robot.y();
                    if(robot.has_orientation())
                        botTheta = robot.orientation();
                    else
                        botTheta = 0;
                    linearTransform(botX, botY, botTheta);
                    if (isTeamYellow) {
                        VisionVelocity::BotPose p1(bsQ.front().first.homeX[botId], bsQ.front().first.homeY[botId], bsQ.front().first.homeTheta[botId]);
                        VisionVelocity::BotPose p2(botX, botY, botTheta);
                        VisionVelocity::calcBotVelocity(p1, p2, timeMs, botVl, botVr);
                        bsMutex->lock();
                        bs->homeVx[botId]= (botX - bs->homeX[botId])/(timeMs*0.001);
                        bs->homeVy[botId]= (botY - bs->homeY[botId])/(timeMs*0.001);
                        bs->homeX[botId] = botX;
                        bs->homeY[botId] = botY;
                        bs->homeTheta[botId] = botTheta;
                        bs->homeVl[botId] = botVl;
                        bs->homeVr[botId] = botVr;
                        bs->homeIsPresent[botId] = true;
                        bsMutex->unlock();
                    } else {
                        VisionVelocity::BotPose p1(bsQ.front().first.awayX[botId], bsQ.front().first.awayY[botId], bsQ.front().first.awayTheta[botId]);
                        VisionVelocity::BotPose p2(botX, botY, botTheta);
                        VisionVelocity::calcBotVelocity(p1, p2, timeMs, botVl, botVr);
                        bsMutex->lock();
                        bs->awayX[botId] = botX;
                        bs->awayY[botId] = botY;
                        bs->awayTheta[botId] = botTheta;
                        bs->awayVl[botId] = botVl;
                        bs->awayVr[botId] = botVr;
                        bs->awayIsPresent[botId] = true;
                        bsMutex->unlock();
                    }
                }                
                bsQ.pop();
                bsMutex->lock();
                bsQ.push(make_pair(*bs, nowTime));
                bsMutex->unlock();
                emit newData();
            }

        }
    }
}


void VisionWorker::onExit()
{
}
