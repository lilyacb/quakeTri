### Earthquake processing and visualization functions

#' plot3compEqData: plot all components of 3-component (Z, N, E) seismic data in one window with labels
#' @param data.df dataframe containing sesimic data for 3 components 
#' @param total.time time interval for the data, default 27s for LEN-DB data
#' @param plotTitle title of the plot
#' @param mag attribute magnitude
#' @param starttime attribute starttime
#' @param sta station code
#' @return ggplot of the seismic data
#' @examples 
#' Usage example
#' @export
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


#' circle_intersection: calculate the intersection of two circles
#' @param center1 a vector containing the center coordinates of the first circle
#' @param radius1 the radius of the first circle
#' @param center2 a vector containing the center coordinates of the second circle
#' @param radius2 the radius of the second circle
#' @return a matrix containing either the ordered pairs of the 2 intersection points or NAs if no intersection is detected
#' @examples 
#' Usage example
#' @export
circle_intersection<-function(center1=c(x1,y1),radius1,center2=c(x2,y2),radius2){
  # Based on the solution of:
  #   x^2 + y^2 + Ax + By + c = 0
  #   x^2 + y^2 + Dx + Ey + F = 0
  #
  # Solution used: Algebraic -- subtract the equations and solve
  # A vector-based solution might also be derived which could be "simpler" than the algebraic derivation
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
    # tangency check? or just return the same point
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


#' circle: calculate a set of points that lie on a circle with a given center and diameter
#' @param center a vector containing the center coordinates of the circle
#' @param diameter the diameter of the circle
#' @param npoints the number of points to generate for the circle, default=100
#' @return a dataframe containing the specified number of points that lie on the circle
#' @examples 
#' Usage example
#' @export 
circle <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#' geo_circle_intersections: use points on two circles (can be created with destPoint) to approximate their intersection points by
#' finding points on each circle that are close together (using distGeo for accurate distance). This functino uses 
#' hclust to separate out the two intersections by geographical distance.
#' @param circle1 dataframe of points on a circle with cols (lon,lat)
#' @param circle2 dataframe of points on a circle with cols (lon,lat)
#' @param intersection_name character string to label the intersection points
#' @return a dataframe of cols (lon,lat) and 2 rows containing the 2 intersection points of the circles (as averages of all 
#' the points in a cluster), if they exist
#' @examples 
#' Usage example
#' RLOK_ELIS_ints<-geo_circle_intersections(RLOK_circle,ELIS_circle,"RLOK_ELIS")
#' @export
geo_circle_intersections<-function(circle1,circle2,intersection_name){
  int1<-c()
  int2<-c()
  ints<-c()
  ints.df<-c()
  for(i in seq(1:dim(circle1)[1])){ #write parallel func???!!!
   # identifier_1_files_ind<-pvec(seq_along(files),function(i)
    #  grepl(identifier_1,files[i],fixed=T),mc.cores=cores)
    for(j in seq(1:dim(circle2)[1])){
      dg<-distGeo(circle1[i,],circle2[j,])/1000<=2.5 # within 2.5 km
      if(dg){
        print("found intersection")
        int<-rbind(circle1[i,],circle2[j,])
        ints<-rbind(ints,int)
      }
    }
  }
  
  # check for contents of ints--if length=0 then exit
  if(length(ints)==0){ 
    stop("No intersections found")
  }
  
  # change ints to a dataframe
  ints.df<-as.data.frame(ints)
  colnames(ints.df)<-c("intersection_longitude","intersection_latitude")
  #print(ints.df)
  
  # parse out the different intersections (some sets of points may represent the same intersection)
  # find where the intersections are different. Create a distance matrix to determine geographical distance
  int_dist<-distm(ints.df,ints.df)/1000 #km
  int_dist.df<-as.data.frame(int_dist)
  # names need to be different
  rnames<-c()
  for(i in seq(1:nrow(int_dist.df))){
    rnames<-c(rnames,paste(intersection_name,"_",i,sep=""))
  }
  #print(rnames)
  rownames(int_dist.df)<-rnames
  colnames(int_dist.df)<-rnames
  int_dist.dist<-as.dist(int_dist.df)
  
  # hclust to find clusters of points
  clust<-hclust(int_dist.dist)
  clust$labels<-rnames
  int_dist.dendro<-clust %>% as.dendrogram() %>% set("labels_to_character")
  plot(int_dist.dendro)
  # grab each clade--get the names and use names to grab the points in the intersection dataframe
  geoCuts<-cutree(int_dist.dendro,h=10)
  # look at the cluster groupings
  print(table(geoCuts))
  # get first group
  first_group<-names(which(geoCuts==1))
  #print(first_group)
  # get second group---add check for tangenecy??
  second_group<-names(which(geoCuts==2))
  #print(second_group)
  # take the average of each set of (lon, lat) (could do something else?)
  # grab the correct points from the int.df that contains the coords
  cutoff<-length(first_group)
  first_group_ints<-ints.df[1:cutoff,]
  print(first_group_ints)
  #print(dim(first_group_ints))
  second_group_ints<-ints.df[(cutoff+1): nrow(ints.df),]
  print(second_group_ints)
  # find the first average
  first_lon_avg<-mean(first_group_ints[,1])
  first_lat_avg<-mean(first_group_ints[,2])
  final_intersection_coords1<-(c(first_lon_avg,first_lat_avg))
  # find the second average
  second_lon_avg<-mean(second_group_ints[,1])
  second_lat_avg<-mean(second_group_ints[,2])
  final_intersection_coords2<-(c(second_lon_avg,second_lat_avg))
  final_int.df<-as.data.frame(rbind(final_intersection_coords1,final_intersection_coords2))
  rownames(final_int.df)<-c(paste(intersection_name,".1",sep=""),paste(intersection_name,".2",sep="")) #change to intersection names
  colnames(final_int.df)<-c("lon","lat")
  print(final_int.df)
  return(final_int.df)
}
