from eventHandling.peak import Peak

class Hit:
    def __init__(self, channel, Peak):
        self.channel = channel
        self.channelNumber = int(self.channel.split("chn")[1])
        self.Peak = Peak
        self.charge = self.Peak.integral_value
        self.hitTime = self.Peak.peak