time1<-1:100
time2<-101:200
time3<-201:300
time4<-301:400
times<-c(time1,time2, time3, time4)
ts1p1<-1/8*sin(2*pi*time1/50)
ts1p2<-0*time2
ts1p3<-0*time3
ts1p4<-0*time4
ts1<-c(ts1p1,ts1p2, ts1p3, ts1p4)
ts<-ts1

plot(ts(ts))

ts2p1<-0*time1
ts2p3<-0*time3
ts2p4<-0*time4
ts2p2<-1/4*sin(2*pi*time2/50)
ts2<-c(ts2p1,ts2p2, ts2p3, ts2p4)
ts<-ts1+ts2

plot(ts(ts))

ts3p1<-0*time1
ts3p2<-0*time2
ts3p4<-0*time4
ts3p3<-1/2*sin(2*pi*time2/50)
ts3<-c(ts3p1,ts3p2, ts3p3, ts3p4)
ts<-ts1+ts2+ts3

plot(ts(ts))

ts4p1<-0*time1
ts4p2<-0*time2
ts4p3<-0*time4
ts4p4<-3/4*sin(2*pi*time2/50)
ts4<-c(ts4p1,ts4p2, ts4p3, ts4p4)
ts<-ts1+ts2+ts3+ts4

plot(ts(ts))




#Then add a sine wave of amplitude 1 and period 8 that operates for t = 101, . . . , 200 but before that is absent.
ts2p1<-0*time1
ts2p2<-sin(2*pi*time2/8)
ts2<-c(ts2p1,ts2p2)
ts<-ts+ts2

plot(ts(ts))
