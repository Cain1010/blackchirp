#include "virtualflowcontroller.h"

#include "virtualinstrument.h"

VirtualFlowController::VirtualFlowController(QObject *parent) : FlowController(parent)
{
    d_subKey = QString("virtual");
    d_prettyName = QString("Virtual Flow Controller");

    d_comm = new VirtualInstrument(d_key,this);
    connect(d_comm,&CommunicationProtocol::logMessage,this,&VirtualFlowController::logMessage);
    connect(d_comm,&CommunicationProtocol::hardwareFailure,[=](){ emit hardwareFailure(); });
}

VirtualFlowController::~VirtualFlowController()
{
}



bool VirtualFlowController::testConnection()
{
    p_readTimer->stop();
    readAll();
    p_readTimer->start();

    emit connected();
    return true;
}

void VirtualFlowController::initialize()
{
    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());
    s.beginGroup(d_key);

    s.beginReadArray(QString("channels"));
    for(int i=0;i<BC_FLOW_NUMCHANNELS;i++)
    {
        s.setArrayIndex(i);
        d_config.add(0.0,s.value(QString("name"),QString("")).toString());
    }
    s.endArray();
    s.endGroup();

    testConnection();
}

Experiment VirtualFlowController::prepareForExperiment(Experiment exp)
{
    return exp;
}

void VirtualFlowController::beginAcquisition()
{
}

void VirtualFlowController::endAcquisition()
{
}

double VirtualFlowController::setFlowSetpoint(const int ch, const double val)
{
    if(ch<0 || ch >= d_config.size())
        return -1.0;

    d_config.set(ch,FlowConfig::Setpoint,val);
    return readFlowSetpoint(ch);
}

double VirtualFlowController::setPressureSetpoint(const double val)
{
    d_config.setPressureSetpoint(val);
    return readPressureSetpoint();
}

double VirtualFlowController::readFlowSetpoint(const int ch)
{
    if(ch < 0 || ch >= d_config.size())
        return -1.0;

    emit flowSetpointUpdate(ch,d_config.setting(ch,FlowConfig::Setpoint).toDouble());
    return d_config.setting(ch,FlowConfig::Setpoint).toDouble();
}

double VirtualFlowController::readPressureSetpoint()
{
    emit pressureSetpointUpdate(d_config.pressureSetpoint());
    return d_config.pressureSetpoint();
}

double VirtualFlowController::readFlow(const int ch)
{
    if(ch < 0 || ch >= d_config.size())
        return -1.0;

    double sp = d_config.setting(ch,FlowConfig::Setpoint).toDouble();
    double noise = sp*((double)(qrand()%100)-50.0)/1000.0;
    double flow = sp + noise;
    d_config.set(ch,FlowConfig::Flow,flow);

    emit flowUpdate(ch,d_config.setting(ch,FlowConfig::Flow).toDouble());
    return d_config.setting(ch,FlowConfig::Setpoint).toDouble();
}

double VirtualFlowController::readPressure()
{
    d_config.setPressure(d_config.pressureSetpoint());

    emit pressureUpdate(d_config.pressure());
    return d_config.pressure();
}

void VirtualFlowController::setPressureControlMode(bool enabled)
{
    d_config.setPressureControlMode(enabled);
    readPressureControlMode();
}

bool VirtualFlowController::readPressureControlMode()
{
    emit pressureControlMode(d_config.pressureControlMode());
    return d_config.pressureControlMode();
}
