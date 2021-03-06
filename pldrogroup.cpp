#include "pldrogroup.h"

PldroGroup::PldroGroup(QObject *parent) :
    Synthesizer(parent)
{
    d_subKey = QString("pldro");
    d_prettyName = QString("PLDRO Oscillators");
    d_commType = CommunicationProtocol::Custom;
    d_threaded = false;

    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());
    s.beginGroup(d_key);
    s.beginGroup(d_subKey);
    //allow hardware limits to be made in settings
    d_minFreq = s.value(QString("minFreq"),0.0).toDouble();
    d_maxFreq = s.value(QString("maxFreq"),1000000.0).toDouble();
    d_txFreq = s.value(QString("lastTxFreq"),0.0).toDouble();
    d_rxFreq = s.value(QString("lastRcvrFreq"),0.0).toDouble();
    //write the settings if they're not there
    s.setValue(QString("lastTxFreq"),d_txFreq);
    s.setValue(QString("lastRcvrFreq"),d_rxFreq);
    s.setValue(QString("minFreq"),d_minFreq);
    s.setValue(QString("maxFreq"),d_maxFreq);
    s.endGroup();
    s.endGroup();
}



bool PldroGroup::testConnection()
{
    readTxFreq();
    readRxFreq();

    emit connected(true);
    return true;
}

void PldroGroup::initialize()
{
    p_comm->initialize();
    testConnection();
}

Experiment PldroGroup::prepareForExperiment(Experiment exp)
{
    return exp;
}

void PldroGroup::beginAcquisition()
{
}

void PldroGroup::endAcquisition()
{
}

void PldroGroup::readTimeData()
{
}

double PldroGroup::readSynthTxFreq()
{
    return d_txFreq;
}

double PldroGroup::readSynthRxFreq()
{
    return d_rxFreq;
}

double PldroGroup::setSynthTxFreq(const double f)
{
    d_txFreq = f;

    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());
    s.beginGroup(d_key);
    s.beginGroup(d_subKey);
    s.setValue(QString("lastTxFreq"),d_txFreq);
    s.endGroup();
    s.endGroup();

    return readTxFreq();
}

double PldroGroup::setSynthRxFreq(const double f)
{
    d_rxFreq = f;

    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());
    s.beginGroup(d_key);
    s.beginGroup(d_subKey);
    s.setValue(QString("lastRcvrFreq"),d_rxFreq);
    s.endGroup();
    s.endGroup();

    return readRxFreq();
}
