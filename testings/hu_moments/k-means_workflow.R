#### K-MEANS FILTERING ###################################################

source("~/Library/Mobile Documents/com~apple~CloudDocs/Analys mikroskopi/k-means_functions.R")

# Define custom classification variables (optional; selected morphology variables used by default)
class.variables = NULL
#class.variables <- c('a.tot','ellipse.perim','perim','maj.axis','min.axis','sphere.vol')

# Set number k-mean classes
k_count <- 10

# Set max number of iterations (*)
max_iter <- 50 
# (*) The function k.means.classification runs the algorithm iteratively,
# and keeps track of how many cells changed class between each iteration.
# If for any iteration there were no class updates, the function considers
# the classification complete, and breaks the run. Otherwise, it continues
# the classification until "max.iter" iterations have been completed.

# Run k-means classification; returns Cell-ID object with column "k",
# indicating the assigned class for each row/cell
X <- k.means.classification(X, k=k_count, max.iterations=max_iter, vars=class.variables)

# Loop over all classes and display "N" random cells
for(i in 1:k_count) show(cimage(X,~cell,channel="BF.out",N=49,subset=k==i,box.size=40))

# Select a class and show histogram of distances to the class centroid
k_select = 9
ggplot(data=subset(X$data,k==k_select & qc==TRUE))+
  geom_histogram(aes(x=k.dist), binwidth=1)

# Set a distance limit, and display only cells more than "k_lim"
# from the class centroid
k_lim = 0
cimage(X,~cell,channel="BF.out",N=49,subset=k==k_select & k.dist>k_lim,box.size=40)

# Filter based on selected class and distance
X<-qc_filter(X,!(k==k_select & k.dist>k.lim))
X<-qc_undo(X)
X<-qc_execute(X)
