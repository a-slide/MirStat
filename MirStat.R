#! /usr/bin/Rscript

options(width = 160) # wider terminal size
library(FactoMineR)
library(RColorBrewer)
library(gplots)
#library(rgl)
x11()

########################################################################################################################
# You need to provide the filename as command line arg
########################################################################################################################

# Inport data in a dataframe
filename = commandArgs(T)[1]
data = read.table(filename, dec=",", sep="\t", header = T, row.names=1)

# Create a list of group names to analyze data separatly
v1 = c("G1-R1", "G1-R2", "G1-R3", "G2-R1", "G2-R2", "G2-R3")
v2 = c("G2-R1", "G2-R2", "G2-R3", "G3-R1", "G3-R2", "G3-R3")
v3 = c("G3-R1", "G3-R2", "G3-R3", "G4-R1", "G4-R2", "G4-R3")
v4 = c("G2-R1", "G2-R2", "G2-R3", "G4-R1", "G4-R2", "G4-R3")
v5 = c("G1-R1", "G1-R2", "G1-R3", "G2-R1", "G2-R2", "G2-R3", "G3-R1", "G3-R2", "G3-R3", "G4-R1", "G4-R2", "G4-R3")

# Define group list and group names
group.list = list(v1,v2,v3,v4,v5)
names(group.list) = c("G1 vs G2", "G2 vs G3", "G3 vs G4", "G2 vs G4", "all samples")

# Analyse each group iteratively
for (group.name in names(group.list)){

    print (paste("ANALYSING GROUP ", group.name))

    #Extract informations for the current group
    group.data = data[group.list[[group.name]],]
    print ("INDIVIDUAL ANALYZED")
    print (group.list[[group.name]])

    # Open a file for writting a report of eliminated and retained mir
    filename = paste("Report_", group.name,".txt", collapse="")

    # Eliminate due to at least 1 NA in the group
    write ("ELIMINATED MIR DUE TO TECHNICAL FILTERS", filename)
    mir_list = paste(names(group.data[,colSums(is.na(group.data)) > 0]), collapse = "\t")
    write (mir_list, filename, append=T)
    group.data = group.data[,colSums(is.na(group.data)) == 0]

    # Eliminate due to all undetermined in the group
    write ("ELIMINATED MIR DUE TO UNDETERMINED VALUES", filename, append=T)
    mir_list = paste(names(group.data[,colSums(group.data !=0) == 0]), collapse = "\t")
    write (mir_list, filename, append=T)
    group.data = group.data[,colSums(group.data !=0) > 0]

    # Retained Mir for analysis
    write ("RETAINED MIR", filename, append=T)
    mir_list = paste(names(group.data), collapse = "\t")
    write (mir_list, filename, append=T)

    # Reduce and center dataset
    group.data = scale(group.data)
    dist_mat = dist(group.data, method = "euclidean")

    # Perform a CAH
    tree = (hclust (dist_mat, method="ward"))
    plot(tree)
    n_clust = nrow(group.data)/3
    rect.hclust (tree, k=n_clust, border=c(1:n_clust))
    dev.print(file = paste("CAH", group.name,".svg", collapse=""), device=svg)

    # Extract classes from the CAH
    class = cutree(tree, k = n_clust)

    # Perform PCA with FactoMineR pakage
    PCA.res = PCA( group.data, scale.unit=T, graph = F)

    #plot3d(PCA.res$ind$coord[,1:3], col=class)
    #text3d(PCA.res$ind$coord[,1:3],texts=rownames(group.data))

    # Create and export plots for all PCA results
    plot(PCA.res, axes = c(1, 2), choix = "ind", title = group.name, col.ind = class)
    dev.print(file = paste("PCA_samples_PC1_PC2", group.name,".svg", collapse=""), device=svg)

    plot(PCA.res, axes = c(1, 3), choix = "ind", title = group.name, col.ind = class)
    dev.print(file = paste("PCA_samples_PC1_PC3", group.name,".svg", collapse=""), device=svg)

    plot(PCA.res, axes = c(1, 2), choix = "var", title = group.name, cex = 0.2, lim.cos2.var = 0.8)
    dev.print(file = paste("PCA_variables_PC1_PC2", group.name,".svg", collapse=""), device=svg)

    plot(PCA.res, axes = c(1, 3), choix = "var", title = group.name, cex = 0.2, lim.cos2.var = 0.8)
    dev.print(file = paste("PCA_variables_PC1_PC3", group.name,".svg", collapse=""), device=svg)

    # Plot the variance explained
    eigen_plot = barplot(PCA.res$eig[,2], names=paste("Dim",1:nrow(PCA.res$eig)),main = group.name)
    dev.print(file = paste("PCA_variance_", group.name,".svg", collapse=""), device=svg)

    # Create summary of highly corelated variables with dim1, dim2, dim3 and export to CSV files
    Dim = dimdesc(PCA.res, axes = c(1, 2, 3))
    output = rbind (
        c("Dim1.Correlation", "Dim1.p-Value"), Dim$Dim.1$quanti,
        c("Dim2.Correlation", "Dim2.p-Value"), Dim$Dim.2$quanti,
        c("Dim3.Correlation", "Dim3.p-Value"), Dim$Dim.3$quanti)
    write.table(output, paste("PCA_var_dim_", group.name,".csv", collapse=""), sep="\t", col.names = F)

    # Write the coordinates of individual in a csv file
    col_names = paste("\t",colnames(PCA.res$ind$coord), collapse="")
    file_out = paste("PCA_ind_coord_", group.name,".csv", collapse="")
    write(col_names, file_out)
    write.table(PCA.res$ind$coord, file_out, sep="\t", append = T, col.names = F, quote = F)

    ## Row clustering
    dist_mat = as.dist(1-cor(group.data, method="pearson"))
    hr <- hclust(dist_mat, method="complete")
    group.data = t(group.data)

    ## Plot heatmap
    heatmap.2(
        group.data,
        Rowv=as.dendrogram(hr),
        Colv= F,
        dendrogram = "row",
        scale="row",
        density.info="density",
        trace = "column",
        tracecol = "black",
        col = colorRampPalette(c("springgreen4","palegreen1","white","coral1","red4")),
        main = group.name,
        cexCol = 0.7,
        cexRow = 0.2)

    dev.print(file = paste("Heatmap_", group.name,".svg", collapse=""), device=svg)
}

warnings()
