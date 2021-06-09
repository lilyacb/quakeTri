### Earthquake processing/visualization functions

# function to plot 3 component earthquake time data (Z,N,E)
plot3compEqData<-function(data.df,total.time=27,plotTitle="plot",mag="mag",starttime="starttime",sta="sta"){ # time in benchmark dataset always 27s
  time.length<-dim(data.df[1:3][])[1]
  time.vec<-as.data.frame(seq(0,total.time,length=time.length))
  colnames(time.vec)<-c("time")
  
  dataZ<-as.data.frame(cbind(time.vec,data.df[,1]))
  colnames(dataZ)<-c("time","Z")
  pZ<-ggplot(data=dataZ,aes(x=time,y=Z))+geom_line()
  
  dataN<-as.data.frame(cbind(time.vec,data.df[,2]))
  colnames(dataN)<-c("time","N")
  pN<-ggplot(data=dataN,aes(x=time,y=N))+geom_line()
  
  
  dataE<-as.data.frame(cbind(time.vec,data.df[,3]))
  colnames(dataE)<-c("time","E")
  pE<-ggplot(data=dataE,aes(x=time,y=E))+geom_line()
  
  grid.arrange(pZ,pN,pE,top=plotTitle,bottom=paste("sta: ", sta, " mag: ", mag, " start: ", starttime," GMT",sep=""))
}

#' cha: convex hull area
#' @param x dataframe of longitudes
#' @param y dataframe of latitudes
#' @return area of the polygon
#' @examples 
#' Usage example
#' @export 
cha<-function(x,y){
  i<-chull(x,y)
  return(areapl(cbind(x[i],y[i])))
}


# circle_intersection, a function to calculate the intersection of two circles given the center coordinates and radii
# Based on the solution of:
#   x^2 + y^2 + Ax + By + c = 0
#   x^2 + y^2 + Dx + Ey + F = 0
#
# Solution used: Algebraic -- subtract the equations and solve
# A vector-based solution may also be derived (it would be shorter, possibly better), but am
#  having issues generalizing the angles needed to define the x and y components 
circle_intersection<-function(center1=c(x1,y1),radius1,center2=c(x2,y2),radius2){
  xInt1=0; xInt2=0; yInt1=0; yInt2=0
  intersection_coords<-c(c(xInt1,yInt1),c(xInt2,yInt2))
  # check if the circles intersect
  x1<-center1[1]
  x2<-center2[1]
  y1<-center1[2]
  y2<-center2[2]
  
  intersection_exists<-((x2-x1)^2 + (y2-y1)^2)^0.5 <= (radius1+radius2)
  if(intersection_exists){
    print("intersection detected")
    # tangency check?
  }
  else{
    print("no intersection detected, aborting calculation")
    return(intersection_coords<-c(c(NA,NA),c(NA,NA)))
  }
  # standard form constants
  A<-(-2*x1)
  B<-(-2*y1)
  C<-(x1^2+y1^2-radius1^2)
  D<-(-2*x2)
  E<-(-2*y2)
  f<-(x2^2+y2^2-radius2^2)
  
  a<-B^2-2*B*E+E^2+(A-D)^2
  b<-B*(-2*f+2*C+(A-D)^2-A^2+A*D)+E*(-2*C+2*f+A^2-A*D)
  c<-f*(f-2*C+A^2-A*D)+C*(C+(A-D)^2-A^2+A*D)
  
  y.a<-(-b+sqrt(b^2-4*a*c))/(2*a)
  y.b<-(-b-sqrt(b^2-4*a*c))/(2*a)
  
  x.a<-(-B*y.a+E*y.a+f-C)/(A-D)
  x.b<-(-B*y.b+E*y.b+f-C)/(A-D)
  
  intersection_coords<-c(c(x.a,y.a),c(x.b,y.b))
  return(intersection_coords)
}

# function to calculate a set of points that lie on a circle with a given center and radius
circle <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}