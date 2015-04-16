#ifndef PULSEGENCONFIG_H
#define PULSEGENCONFIG_H

#include <QSharedDataPointer>
#include <QList>
#include <QVariant>
#include <QMap>

class PulseGenConfigData;

class PulseGenConfig
{
public:
    PulseGenConfig();
    PulseGenConfig(const PulseGenConfig &);
    PulseGenConfig &operator=(const PulseGenConfig &);
    ~PulseGenConfig();

    enum ActiveLevel { ActiveLow, ActiveHigh };
    enum Setting { Delay, Width, Enabled, Level, Name };

    struct ChannelConfig {
        int channel;
        QString channelName;
        bool enabled;
        double delay;
        double width;
        ActiveLevel level;

        ChannelConfig() : channel(-1), enabled(false), delay(-1.0), width(-1.0), level(ActiveHigh) {}
    };

    PulseGenConfig::ChannelConfig at(const int i) const;
    int size() const;
    bool isEmpty() const;
    QVariant setting(const int index, const PulseGenConfig::Setting s) const;
    PulseGenConfig::ChannelConfig settings(const int index) const;
    QMap<QString,QPair<QVariant,QString>> headerMap() const;

    void set(const int index, const PulseGenConfig::Setting s, const QVariant val);
    void set(const int index, const PulseGenConfig::ChannelConfig cc);
    void add(const QString name, const bool enabled, const double delay, const double width, const ActiveLevel level);

private:
    QSharedDataPointer<PulseGenConfigData> data;
};

class PulseGenConfigData : public QSharedData
{
public:
    QList<PulseGenConfig::ChannelConfig> config;
};

Q_DECLARE_METATYPE(PulseGenConfig::ActiveLevel)

#endif // PULSEGENCONFIG_H