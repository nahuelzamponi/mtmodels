
library(igraph)

#########################functions#########################
#mod
#b = a - m.*floor(a./m)
mod <- function(x,m) {
	t1<-floor(x/m)
	return(x-t1*m)
}

#network operators
ci <- function(g) {
	pick <- sample(N, 1) #pick a node
	Prob <- round(runif(2), 2)

	nlink1 <- nearlink1[pick]
	nlink2 <- nearlink2[pick]
	flink <- farlink[pick]

	if (p1 < Prob[1]) { g[pick, nlink1] <- 0
		g[pick, nlink2] <- 0 }
	else { g[pick, nlink1] <- 1
		g[pick, nlink2] <- 1 }

	if (p2 < Prob[2]) { g[pick, flink] <- 0 }
	else { g[pick, flink] <- 1 }

	return(g) }

#########################parameters#########################
#size of the array
L = 20 #please, even number...
N = L^2 #number of nodes

T = 1000 #iterations

p1 <- 0.65
p2 <- 0.5

#node labeling, 1 to N
x <- vector()
y <- vector()
for (i in 1:N) {
	x[i] <- mod(i, L)
	if (x[i] == 0) { x[i] <- L }
	y[i] <- ceiling(i/L)
}

#neighbors
#tip-to-tip
nearlink1<-rep(0, N)
nearlink2<-rep(0, N)
for (i in 1:N) {
	if (i <= N) {
		nearlink1[i] <- i+1
	}
	if (i > 1) {
		nearlink2[i] <- i-1
	}
}
nearlink1[N]=1
nearlink2[1]=N

#tip-to-side
#it generates a vector in which if farlink[i]=0 then the condition farlink[j]=i holds
farlink <- rep(0, N)
for (i in 1:(N-L)) {
	if (((y[i]%%2 != 0) == 'TRUE') && ((x[i]%%2 != 0) == 'TRUE')) { #odd horizontal rows -> bonds 1, 3, 5, ...
		farlink[i] <- i+L
	}
	if (((y[i]%%2 != 0) == 'FALSE') && ((x[i]%%2 != 0) == 'FALSE')) { #even horizontal rows -> bonds 2, 4, 6, ...
		farlink[i] <- i+L
	}
}

if ((L %% 2 != 0) == 'TRUE') { pos <- seq((N-L+1), N, 2) #odd
	print('odd') }
if ((L %% 2 != 0) == 'FALSE') { pos <- seq((N-L+2), N, 2) #even 
	print('even') }
cont=2
for (i in 1:length(pos)) { farlink[[pos[i]]] <- cont
				cont=cont+2 }
aux <- c()
for (i in 1:length(farlink)) {
	if (farlink[i] == 0) { aux <- c(aux, which(farlink == i)) }
	else { aux <- c(aux, farlink[i]) }
}
farlink <- aux

########################################################################################################################
######################################## simulation ####################################################################
########################################################################################################################

Ng <- c()
N2 <- c()
ms <- c()
S <- c()
NgN <- c()
mk <- c()

#empty network
g <- graph.empty(n=N, directed=FALSE)

for (t in 1:T) {
	g <- ci(g)
	g <- simplify(g, remove.multiple=T, remove.loops=T)
	clg <- clusters(g)

	sizes <- sort(clg$csize, decreasing=TRUE)

	Ng <- c(Ng, sizes[1])
	N2 <- c(N2, sizes[2])
	
	#ms
	sum_num <- c()
	sum_den <- c()
	s <- sizes[2:length(sizes)] # cluster sizes excluding Ng
	us <- unique(s) # unique cluster sizes (bins)
	for (kk in 1:length(us)) {
		ns <- length(s[which(s == us[kk])])/length(s) # proportion of clusters of size si
		sum_num <- c(sum_num, (us[kk]^2)*ns) # si^2 x ns
		sum_den <- c(sum_den, us[kk]*ns) # si x ns
	}
	ms <- c(ms, sum(sum_num)/sum(sum_den))

	S <- c(S, length(sizes))
	NgN <- c(NgN, sizes[1]/sum(sizes))
	mk <- c(mk, mean(degree(g)))

	print(paste('iterations =', t))
}



#plot quantities
par(mfrow=c(2, 3))
plot(Ng ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='Ng')
plot(N2 ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='N2')
plot(ms ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='<s>')
plot(S ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='S')
plot(NgN ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='Ng/N')
plot(mk ~ c(1:T), pch=21, bg='white', type='p', las=1, xlab='iterations', ylab='', main='<k>')


#plot graph
V(g)$label <- NA
V(g)$size <- degree(g)*.5
E(g)$width <- 2
#plot(g, edge.color="firebrick1", vertex.color="green", layout=layout_with_kk)
