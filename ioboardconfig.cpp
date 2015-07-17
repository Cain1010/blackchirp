#include "ioboardconfig.h"

#include <QStringList>
#include <QSettings>
#include <QApplication>

class IOBoardConfigData : public QSharedData
{
public:
    QMap<int,QPair<bool,QString>> analog;
    QMap<int,QPair<bool,QString>> digital;

    int numAnalog;
    int numDigital;
    int reservedAnalog;
    int reservedDigital;

};

IOBoardConfig::IOBoardConfig() : data(new IOBoardConfigData)
{
    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());
    s.beginGroup(QString("ioboard"));
    s.beginGroup(s.value(QString("subKey"),QString("virtual")).toString());

    data->numAnalog = qBound(0,s.value(QString("numAnalog"),4).toInt(),16);
    data->numDigital = qBound(0,s.value(QString("numDigital"),16-data->numAnalog).toInt(),16);
    data->reservedAnalog = qMin(data->numAnalog,s.value(QString("reservedAnalog"),0).toInt());
    data->reservedDigital = qMin(data->numDigital,s.value(QString("reservedDigital"),0).toInt());

    s.endGroup();
    s.endGroup();

    s.beginGroup(QString("iobconfig"));
    s.beginReadArray(QString("analog"));
    for(int i=0; i<data->numAnalog-data->reservedAnalog; i++)
    {
        s.setArrayIndex(i);
        QString name = s.value(QString("name"),QString("")).toString();
        bool enabled = s.value(QString("enabled"),false).toBool();
        data->analog.insert(i,qMakePair(enabled,name));
    }
    s.endArray();
    s.beginReadArray(QString("digital"));
    for(int i=0; i<data->numDigital-data->reservedDigital; i++)
    {
        s.setArrayIndex(i);
        QString name = s.value(QString("name"),QString("")).toString();
        bool enabled = s.value(QString("enabled"),false).toBool();
        data->digital.insert(i,qMakePair(enabled,name));
    }
    s.endArray();

    s.endGroup();
}

IOBoardConfig::IOBoardConfig(const IOBoardConfig &rhs) : data(rhs.data)
{

}

IOBoardConfig &IOBoardConfig::operator=(const IOBoardConfig &rhs)
{
    if (this != &rhs)
        data.operator=(rhs.data);
    return *this;
}

IOBoardConfig::~IOBoardConfig()
{

}

void IOBoardConfig::setAnalogChannel(int ch, bool enabled, QString name)
{
    if(data->analog.contains(ch))
        data->analog[ch] = qMakePair(enabled,name);
    else
        data->analog.insert(ch,qMakePair(enabled,name));
}

void IOBoardConfig::setDigitalChannel(int ch, bool enabled, QString name)
{
    if(data->digital.contains(ch))
        data->digital[ch] = qMakePair(enabled,name);
    else
        data->digital.insert(ch,qMakePair(enabled,name));
}

void IOBoardConfig::setAnalogChannels(const QMap<int, QPair<bool, QString> > l)
{
    data->analog = l;
}

void IOBoardConfig::setDigitalChannels(const QMap<int, QPair<bool, QString> > l)
{
    data->digital = l;
}

int IOBoardConfig::numAnalogChannels() const
{
    return data->numAnalog;
}

int IOBoardConfig::numDigitalChannels() const
{
    return data->numDigital;
}

int IOBoardConfig::reservedAnalogChannels() const
{
    return data->reservedAnalog;
}

int IOBoardConfig::reservedDigitalChannels() const
{
    return data->numDigital;
}

bool IOBoardConfig::isAnalogChEnabled(int ch) const
{
    if(data->analog.contains(ch))
        return data->analog.value(ch).first;
    else
        return false;
}

bool IOBoardConfig::isDigitalChEnabled(int ch) const
{
    if(data->digital.contains(ch))
        return data->digital.value(ch).first;
    else
        return false;
}

QMap<int, QPair<bool, QString> > IOBoardConfig::analogList() const
{
    return data->analog;
}

QMap<int, QPair<bool, QString> > IOBoardConfig::digitalList() const
{
    return data->digital;
}

QMap<QString, QPair<QVariant, QString> > IOBoardConfig::headerMap() const
{
    QMap<QString, QPair<QVariant, QString> > out;

    auto it = data->analog.constBegin();
    QString prefix = QString("IOBoardConfig");
    QString empty = QString("");
    for(;it != data->analog.constEnd(); it++)
    {
        out.insert(prefix+QString("Analog.")+QString::number(it.key()+data->reservedAnalog)+QString(".Enabled"),qMakePair(it.value().first,empty));
        out.insert(prefix+QString("Analog.")+QString::number(it.key()+data->reservedAnalog)+QString(".Name"),qMakePair(it.value().second,empty));
    }
    it = data->digital.constBegin();
    for(;it != data->digital.constEnd(); it++)
    {
        out.insert(prefix+QString("Digital.")+QString::number(it.key()+data->reservedDigital)+QString(".Enabled"),qMakePair(it.value().first,empty));
        out.insert(prefix+QString("Digital.")+QString::number(it.key()+data->reservedDigital)+QString(".Name"),qMakePair(it.value().second,empty));
    }

    return out;
}

void IOBoardConfig::saveToSettings() const
{
    QSettings s(QSettings::SystemScope,QApplication::organizationName(),QApplication::applicationName());

    s.beginGroup(QString("iobconfig"));
    s.remove(QString("analog"));
    s.beginWriteArray(QString("analog"));
    for(int i=0; i<data->numAnalog-data->reservedAnalog; i++)
    {
        s.setArrayIndex(i);
        QString name = data->analog.value(i).second;
        if(name.isEmpty())
            name = QString("ain%1").arg(i+data->reservedAnalog);
        s.setValue(QString("name"),name);
        s.setValue(QString("enabled"),data->analog.value(i).first);
    }
    s.endArray();
    s.remove(QString("digital"));
    s.beginWriteArray(QString("digital"));
    for(int i=0; i<data->numDigital-data->reservedDigital; i++)
    {
        s.setArrayIndex(i);
        QString name = data->digital.value(i).second;
        if(name.isEmpty())
            name = QString("din%1").arg(i+data->reservedDigital);
        s.setValue(QString("name"),name);
        s.setValue(QString("enabled"),data->digital.value(i).first);
    }
    s.endArray();

    s.endGroup();
    s.sync();
}
