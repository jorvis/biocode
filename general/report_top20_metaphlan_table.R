#!/usr/bin/env Rscript

# Written by Yang Song

args <- commandArgs( TRUE )

#####install packages-----------------------------
#source("https://bioconductor.org/biocLite.R")

#biocLite('phyloseq')


#install.packages("tidyr",repos='http://cran.rstudio.com/')
#install.packages("ggplot2",repos='http://cran.rstudio.com/')



## ------------------------------------------------------------------------
library(phyloseq);library(tidyr)




#read metaphylan table 
abundance=read.table(args[1],sep="\t",header=T,row.names=1,quote="")


#create the metadata with phenotypes
metadata=data.frame(names(abundance))

#split the names of abundance table name into  phenotype table
foo <- data.frame(do.call('rbind', strsplit(as.character(metadata$names.abundance.),'.',fixed=TRUE)))

#create the type column
metadata$type=paste(foo$X3,foo$X4,sep="_")

#create sample name column
metadata$ID=foo$X1
#names(abundance)=foo$X1

# create the clade dataframe 

tax=data.frame(row.names(abundance))

#rename the column name of the clade dataframe
names(tax)[1]="taxa"

#split the names into 8 taxanomic levels
tax1=separate(tax, taxa, paste0("X",1:8), sep="\\|")
tax1[is.na(tax1)] <- " "


#preapre the 3 required files  for phyloseq objects, OTU table, taxa table and phenotype file

## OTU table
OTU = otu_table(abundance, taxa_are_rows = T)

#tax table
TAX = tax_table(tax1)

rownames(TAX)=rownames(OTU)

#phenotype table
sampledata = sample_data(metadata)

rownames(sampledata)<-colnames(OTU)

#merge 3 datasets into one R object
physeq1=phyloseq(OTU,TAX,sampledata)

#rename the column names of tax_table
colnames(tax_table(physeq1)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","species","strain")


## ------------------------------------------------------------------------

#function to extract different clads from the whole object

subset_taxa=function (physeq, ...) 
{
    if (is.null(tax_table(physeq))) {
        cat("Nothing subset. No taxonomyTable in physeq.\n")
        return(physeq)
    }
    else {
        oldMA <- as(tax_table(physeq), "matrix")
        oldDF <- data.frame(oldMA)
        newDF <- subset(oldDF, ...)
        newMA <- as(newDF, "matrix")
        if (inherits(physeq, "taxonomyTable")) {
            return(tax_table(newMA))
        }
        else {
            tax_table(physeq) <- tax_table(newMA)
            return(physeq)
        }
    }
}

#split the object into different  clade objects

kindom=subset_taxa( physeq1,Phylum==" ")
phylum=subset_taxa(physeq1,Class==" "&Kingdom!=" "&Phylum!=" ")
class=subset_taxa(physeq1,Order==" "&Kingdom!=" "&Phylum!=" "&Class!=" ")
order=subset_taxa(physeq1,Family==" "&Kingdom!=" "&Phylum!=" "&Class!=" "&Order!=" ")
family=subset_taxa(physeq1,Genus==" "&Kingdom!=" "&Phylum!=" "&Class!=" "&Order!=" "&Family!=" ")
genus=subset_taxa(physeq1,species==" "&Kingdom!=" "&Phylum!=" "&Class!=" "&Order!=" "&Family!=" "&Genus!=" ")
species=subset_taxa(physeq1,strain==" "&Kingdom!=" "&Phylum!=" "&Class!=" "&Order!=" "&Family!=" "&Genus!=" "&species!=" ")





## ------------------------------------------------------------------------

#plot richness plot using genus abundance based on genus relative abundance
#pagesToOutputOn = 1     
 #pdf("arg[2]")    
#plot_richness(genus, measures=c("Shannon"))
#p = plot_richness(genus, x = "type")
#print(p + geom_boxplot(data = p$data, aes(x = variable, y = value, color = NULL), alpha = 0.1))

#dev.off()



## ------------------------------------------------------------------------

#plot top20 species
library(ggplot2)
top20 <- names(sort(taxa_sums(species), decreasing=TRUE))[1:20]
ent20 = prune_taxa(top20, species)
pdf("top20.species.pdf")    
#plot the top20 species abunance in one plot
plot_bar(ent20,fill = "species")

#plot the top20 species abundance split by samples and type

p = plot_bar(ent20, "Genus", fill="Genus", facet_grid=type~ID)
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")


dev.off()




## ------------------------------------------------------------------------

#MDS plot based on the species comoposition

#pdf("~/Desktop/cordinate.pdf")    
#plot_ordination(species, ordinate(species, "MDS"),color="ID") + geom_point(size = 5)

#dev.off()



## ------------------------------------------------------------------------

#SAVE DATA

#write.table(data.frame(otu_table(genus)),file="~/Desktop/genus.txt",sep="\t",quote=F)

