% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/earthquakeFunctions.R
\name{geo_circle_intersections}
\alias{geo_circle_intersections}
\title{geo_circle_intersections: use points on two circles (can be created with destPoint) to approximate their intersection points by
finding points on each circle that are close together (using distGeo for accurate distance). This functino uses
hclust to separate out the two intersections by geographical distance.}
\usage{
geo_circle_intersections(circle1, circle2, intersection_name)
}
\arguments{
\item{circle1}{dataframe of points on a circle with cols (lon,lat)}

\item{circle2}{dataframe of points on a circle with cols (lon,lat)}

\item{intersection_name}{character string to label the intersection points}
}
\value{
a dataframe of cols (lon,lat) and 2 rows containing the 2 intersection points of the circles (as averages of all
the points in a cluster), if they exist
}
\description{
geo_circle_intersections: use points on two circles (can be created with destPoint) to approximate their intersection points by
finding points on each circle that are close together (using distGeo for accurate distance). This functino uses
hclust to separate out the two intersections by geographical distance.
}
\examples{
Usage example
RLOK_ELIS_ints<-geo_circle_intersections(RLOK_circle,ELIS_circle,"RLOK_ELIS")
}
