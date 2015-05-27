library(seqLogo)
list <- list.files("C:/Users/David/Desktop/")
counter = 0
for(i in 1:length(list)){
	x<-regexpr("PWM..txt",list[i])
	if(x[1]>=0){
			counter = counter +1
			y<-paste("C:/Users/David/Desktop/",list[i], sep = "")
			m<-read.table(y)
			pwm <- makePWM(m)
			seqLogo(pwm, ic.scale=FALSE)
			
			out<- paste("C:/Users/David/Desktop/PWM",counter,sep = "")
			out<- paste(out,".jpg",sep = "")
			
			savePlot(filename=out, type="jpeg")

		}
}
