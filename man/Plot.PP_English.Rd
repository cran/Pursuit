\name{Plot.PP}

\alias{Plot.PP}

\title{Graphics of the Projection Pursuit (PP).}

\description{Graphics of the Projection Pursuit (PP).}

\usage{Plot.PP(PP, titles = NA, xlabel = NA, ylabel = NA, posleg = 2, boxleg = TRUE,
        size = 1.1, grid = TRUE, color = TRUE, classcolor = NA, linlab = NA, 
        axesvar = TRUE, axes = TRUE, savptc = FALSE, width = 3236, height = 2000, 
        res = 300, casc = TRUE)}
 
\arguments{
  \item{PP}{Data of the PP_Optimizer function.}
  \item{titles}{Titles of the graphics, if not set, assumes the default text.}
  \item{xlabel}{Names the X axis, if not set, assumes the default text.}
  \item{ylabel}{Names the Y axis, if not set, assumes the default text.}
  \item{posleg}{0 with no caption,\cr
                1 for caption in the left upper corner,\cr
                2 for caption in the right upper corner (default),\cr
                3 for caption in the right lower corner,\cr
                4 for caption in the left lower corner.}
  \item{boxleg}{Puts the frame in the caption (default = TRUE).}
  \item{size}{Size of the points in the graphs.}
  \item{grid}{Put grid on graphs (default = TRUE).}
  \item{color}{Colored graphics (default = TRUE).}
  \item{classcolor}{Vector with the colors of the classes.}
  \item{linlab}{Vector with the labels for the observations.}
  \item{axesvar}{Puts axes of rotation of the variables, only when dimproj > 1 (default = TRUE).}
  \item{axes}{Plots the X and Y axes (default = TRUE).}
  \item{savptc}{Saves graphics images to files (default = FALSE).} 
  \item{width}{Graphics images width when savptc = TRUE (defaul = 3236).}
  \item{height}{Graphics images height when savptc = TRUE (default = 2000).}
  \item{res}{Nominal resolution in ppi of the graphics images when savptc = TRUE (default = 300).}
  \item{casc}{Cascade effect in the presentation of the graphics (default = TRUE).}
}

\value{Graph of the evolution of the indices, and graphs whose data were reduced in two dimensions.}

\author{
Paulo Cesar Ossani

Marcelo Angelo Cirillo
}

\seealso{\code{\link{PP_Optimizer}} and \code{\link{PP_Index}}}

\examples{
data(iris) # dataset

# Example 1 - Without the classes in the data
data <- iris[,1:4]

findex <- "kurtosismax" # index function

dim <- 1 # dimension of data projection

sphere <- TRUE # spherical data

res <- PP_Optimizer(data = data, class = NA, findex = findex,
                    optmethod = "GTSA", dimproj = dim, sphere = sphere, 
                    weight = TRUE, lambda = 0.1, r = 1, cooling = 0.9, 
                    eps = 1e-3, maxiter = 500, half = 30)

Plot.PP(res, titles = NA, posleg = 1, boxleg = FALSE, color = TRUE,
        linlab = NA, axesvar = TRUE, axes = TRUE, savptc = FALSE, 
        width = 3236, height = 2000, res = 300, casc = FALSE)


# Example 2 - With the classes in the data
class <- iris[,5] # data class

res <- PP_Optimizer(data = data, class = class, findex = findex,
                    optmethod = "GTSA", dimproj = dim, sphere = sphere, 
                    weight = TRUE, lambda = 0.1, r = 1, cooling = 0.9, 
                    eps = 1e-3, maxiter = 500, half = 30)

tit <- c(NA,"Graph example") # titles for the graphics

Plot.PP(res, titles = tit, posleg = 1, boxleg = FALSE, color = TRUE, 
        classcolor = c("blue3","red","goldenrod3"), linlab = NA, 
        axesvar = TRUE, axes = TRUE, savptc = FALSE, width = 3236,
        height = 2000, res = 300, casc = FALSE)


# Example 3 - Without the classes in the data, but informing 
#             the classes in the plot function
res <- PP_Optimizer(data = data, class = NA, findex = "Moment",
                    optmethod = "GTSA", dimproj = 2, sphere = sphere, 
                    weight = TRUE, lambda = 0.1, r = 1, cooling = 0.9, 
                    eps = 1e-3, maxiter = 500, half = 30)

lin <- c(rep("a",50),rep("b",50),rep("c",50)) # data class

Plot.PP(res, titles = NA, posleg = 1, boxleg = FALSE, color = TRUE,
        linlab = lin, axesvar = TRUE, axes = TRUE, savptc = FALSE, 
        width = 3236, height = 2000, res = 300, casc = FALSE)


# Example 4 - With the classes in the data, but not informed in plot function
class <- iris[,5] # data class

dim <- 2 # dimension of data projection

findex <- "lda" # index function

res <- PP_Optimizer(data = data, class = class, findex = findex,
                    optmethod = "GTSA", dimproj = dim, sphere = sphere, 
                    weight = TRUE, lambda = 0.1, r = 1, cooling = 0.9, 
                    eps = 1e-3, maxiter = 500, half = 30)

tit <- c("",NA) # titles for the graphics

Plot.PP(res, titles = tit, posleg = 1, boxleg = FALSE, color = TRUE,
        linlab = NA, axesvar = TRUE, axes = TRUE, savptc = FALSE, 
        width = 3236, height = 2000, res = 300, casc = FALSE)

}
