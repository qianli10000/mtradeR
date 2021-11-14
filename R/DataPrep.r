
#' @title Data preprocessing
#' @description This function transforms the raw counts OTU table into compositional abundance table based on the specified taxonomy level(s).
#' @param counts A data frame of OTU-level raw counts with rows as taxa and columns as samples. The row names of counts must be the full annotation of each OTU.
#' @param level One or multiple levels in 'kingdom','phylum','class','order','family','genus','species'. Default to 'species'.
#' @param covdata A data frame of subject ID, set ID, sample ID, disease outcome, and covariates to be included in trajectory analysis.
#' @param sample_id Specify the variable for microbiome sample ID in covdata, which should be identical to the column names of counts.
#' @param composition If composition = TRUE, then returns composition. Otherwise, return raw counts in output. Default to TRUE.
#' @param cleananno If cleananno=TRUE, then it cleans out special characters from OTU annotations. Default to FALSE.
#' @return  \item{$relabun}{Relative abundance table with rows as taxa and columns as samples}
#'          \item{$meta_data}{A data frame of all covariates (including disease outcome) to be used in JointMatch or JointSingle}
#' @export


DataPrep <- function(counts,level,covdata,sample_id,composition=TRUE,cleananno=FALSE) {

  if(missing(covdata)==TRUE )
  {stop("Please provide correct meta data")}

  rownames(covdata) <- covdata[,sample_id]

  if(length(intersect(covdata[,sample_id],colnames(counts))) < length(covdata[,sample_id]))
    {stop("Sample IDs dont match with abundance data Sample IDs")}


  if (missing(level) == TRUE || length(level) == 1)

  {


    ## Replace the space after the semicolon
    rownames(counts) <- gsub(pattern = '; ',replacement = ';',rownames(counts))
    ## To check if at least 5 semicolons are present, if not throw an error
    chk <- stringr::str_count(rownames(counts),";")
    if (length(chk[chk<5]) > 0)
      stop("Please ensure data is provided at either genus or species level")


    if(length(chk[chk==6])==length(chk))

    {

      print("The lowest rank is species")

      if(missing(composition)){composition=FALSE}

      #

      # Clean OTU annotation and convert to composition abundance
      levels=as.data.frame(do.call('rbind',strsplit(rownames(counts),';')))

      if(missing(cleananno)){cleananno=FALSE}

      if(cleananno==TRUE)
      {
        levels$V7=gsub(pattern = '\\[',replacement = '',levels$V7)
        levels$V7=gsub(pattern = '\\]',replacement = '',levels$V7)

        levels=as.data.frame(apply(levels,2,function(x){
          y=gsub(' ','.',x)
          y=gsub('-','.',y)

          final_obj <- list(relabun=y,meta_data=covdata)
          return(final_obj)
        }
        )) }

      colnames(levels)=c('kingdom','phylum','class','order','family','genus','species')
      otu.names=paste('k__',levels[,1],'.p__',levels[,2],'.c__',levels[,3],'.o__',levels[,4],'.f__',levels[,5],'.g__',levels[,6],'.s__',levels[,7],sep = '')
      rownames(counts)=otu.names


      if(composition==TRUE)
      {
        counts=t(t(counts)/colSums(counts))
      }

      if(missing(level)){
        final_obj <- list(relabun=counts,meta_data=covdata)
        return(final_obj)}  ## Returns at the species level by default




      # Microbiome composition at each rank

      if (level=="species")
      {
        final_obj <- list(relabun=counts,meta_data=covdata)
        return(final_obj)
      }

      if(level=="genus")

      {
        lv_g=stats::aggregate(counts,by=list(levels$genus),FUN=sum)
        order.id=order(levels$genus)
        rownames(lv_g)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],'.f__',levels[order.id,5],'.g__',levels[order.id,6],sep = ''))

        lv_g <- lv_g[-1]
        final_obj <- list(relabun=lv_g,meta_data=covdata)
        return(final_obj)

      }


      if(level=="family")

      {
        lv_f=stats::aggregate(counts,by=list(levels$family),FUN=sum)

        order.id=order(levels$family)
        rownames(lv_f)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],'.f__',levels[order.id,5],sep = ''))

        lv_f <- lv_f[-1]
        final_obj <- list(relabun=lv_f,meta_data=covdata)
        return(final_obj)

      }

      if(level=="order")

      {
        lv_o=stats::aggregate(counts,by=list(levels$order),FUN=sum)
        order.id=order(levels$order)
        rownames(lv_o)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],sep = ''))


        lv_o <- lv_o[-1]
        final_obj <- list(relabun=lv_o,meta_data=covdata)
        return(final_obj)

      }

      if(level=="class")

      {
        lv_c=stats::aggregate(counts,by=list(levels$class),FUN=sum)
        order.id=order(levels$class)
        rownames(lv_c)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],sep = ''))


        lv_c <- lv_c[-1]
        final_obj <- list(relabun=lv_c,meta_data=covdata)
        return(final_obj)

      }

      if(level=="phylum")

      {
        lv_p=stats::aggregate(counts,by=list(levels$phylum),FUN=sum)
        order.id=order(levels$phylum)
        rownames(lv_p)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],sep = ''))


        lv_p <- lv_p[-1]
        final_obj <- list(relabun=lv_p,meta_data=covdata)
        return(final_obj)

      }
      stop()

    }


    ## If input is only up until genus level

    if(length(chk[chk==5])==length(chk))

    {
      print("The lowest rank is genus")

      if(missing(composition)){composition=FALSE}

      #

      # Clean OTU annotation and convert to composition abundance
      levels=as.data.frame(do.call('rbind',strsplit(rownames(counts),';')))

      if(missing(cleananno)){composition=FALSE}

      if(cleananno==TRUE)
      {
        levels$V6=gsub(pattern = '\\[',replacement = '',levels$V6)
        levels$V6=gsub(pattern = '\\]',replacement = '',levels$V6)
        levels=as.data.frame(apply(levels,2,function(x){
          y=gsub(' ','.',x)
          y=gsub('-','.',y)


          final_obj <- list(relabun=y,meta_data=covdata)
          return(final_obj)
        }
        ))}

      colnames(levels)=c('kingdom','phylum','class','order','family','genus')
      otu.names=paste('k__',levels[,1],'.p__',levels[,2],'.c__',levels[,3],'.o__',levels[,4],'.f__',levels[,5],'.g__',levels[,6],sep = '')
      rownames(counts)=otu.names


      if(composition==TRUE)
      {
        counts=t(t(counts)/colSums(counts))
      }

      if(missing(level)){
        final_obj <- list(relabun=counts,meta_data=covdata)
        return(final_obj)}  ## Returns at the Genus level by default


      # Microbiome composition at each rank


      if(level=="genus")

      {
        lv_g=stats::aggregate(counts,by=list(levels$genus),FUN=sum)
        order.id=order(levels$genus)
        rownames(lv_g)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],'.f__',levels[order.id,5],'.g__',levels[order.id,6],sep = ''))

        lv_g <- lv_g[-1]
        final_obj <- list(relabun=lv_g,meta_data=covdata)
        return(final_obj)

      }


      if(level=="family")

      {
        lv_f=stats::aggregate(counts,by=list(levels$family),FUN=sum)

        order.id=order(levels$family)
        rownames(lv_f)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],'.f__',levels[order.id,5],sep = ''))
        lv_f <- lv_f[-1]
        final_obj <- list(relabun=lv_f,meta_data=covdata)
        return(final_obj)

      }

      if(level=="order")

      {
        lv_o=stats::aggregate(counts,by=list(levels$order),FUN=sum)
        order.id=order(levels$order)
        rownames(lv_o)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                    '.o__',levels[order.id,4],sep = ''))
        lv_o <- lv_o[-1]
        final_obj <- list(relabun=lv_o,meta_data=covdata)
        return(final_obj)

      }

      if(level=="class")

      {
        lv_c=stats::aggregate(counts,by=list(levels$class),FUN=sum)
        order.id=order(levels$class)
        rownames(lv_c)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],sep = ''))

        lv_c <- lv_c[-1]
        final_obj <- list(relabun=lv_c,meta_data=covdata)
        return(final_obj)

      }

      if(level=="phylum")

      {
        lv_p=stats::aggregate(counts,by=list(levels$phylum),FUN=sum)
        order.id=order(levels$phylum)
        rownames(lv_p)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],sep = ''))

        lv_p <- lv_p[-1]
        final_obj <- list(relabun=lv_p,meta_data=covdata)
        return(final_obj)

      }
      stop()

    }

  }


  ## Following code is for rbinding more than one level


  rownames(counts) <- gsub(pattern = '; ',replacement = ';',rownames(counts))
  if(missing(composition)){composition=FALSE}

  ## To check if at least 5 semicolons are present, if not throw an error
  chk <- stringr::str_count(rownames(counts),";")
  if (length(chk[chk<5]) > 0)
    stop("Please ensure data is provided at either genus or species level")

  if(length(chk[chk==6])==length(chk))

  {

    print("Lowest level is species")

    # Clean OTU annotation and convert to composition abundance
    levels=as.data.frame(do.call('rbind',strsplit(rownames(counts),';')))

    if(missing(cleananno)){composition=FALSE}

    if(cleananno==TRUE)
    {
      levels$V7=gsub(pattern = '\\[',replacement = '',levels$V7)
      levels$V7=gsub(pattern = '\\]',replacement = '',levels$V7)
      levels=as.data.frame(apply(levels,2,function(x){
        y=gsub(' ','.',x)
        y=gsub('-','.',y)

        final_obj <- list(relabun=y,meta_data=covdata)
        return(final_obj)
      }
      ))}

    colnames(levels)=c('kingdom','phylum','class','order','family','genus','species')
    otu.names=paste('k__',levels[,1],'.p__',levels[,2],'.c__',levels[,3],'.o__',levels[,4],'.f__',levels[,5],'.g__',levels[,6],'.s__',levels[,7],sep = '')
    rownames(counts)=otu.names


    if(composition==TRUE)
    {
      counts=t(t(counts)/colSums(counts))
    }



    ## Phylum

    lv_p=stats::aggregate(counts,by=list(levels$phylum),FUN=sum)
    order.id=order(levels$phylum)
    rownames(lv_p)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],sep = ''))

    ## Class
    lv_c=stats::aggregate(counts,by=list(levels$class),FUN=sum)
    order.id=order(levels$class)
    rownames(lv_c)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],sep = ''))

    ## Order
    lv_o=stats::aggregate(counts,by=list(levels$order),FUN=sum)
    order.id=order(levels$order)
    rownames(lv_o)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],sep = ''))

    ## Family
    lv_f=stats::aggregate(counts,by=list(levels$family),FUN=sum)
    order.id=order(levels$family)
    rownames(lv_f)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],'.f__',levels[order.id,5],sep = ''))
    ##Genus
    lv_g=stats::aggregate(counts,by=list(levels$genus),FUN=sum)
    order.id=order(levels$genus)
    rownames(lv_g)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],'.f__',levels[order.id,5],'.g__',levels[order.id,6],sep = ''))


    joined_counts <- NULL#matrix(nrow = 1,ncol = dim(lv_g)[2],dimnames = list(c(),c(colnames(lv_g))))

    for(i in 1:length(level))
    {

      if (level[i]=="species")
      {
        joined_counts <- rbind(joined_counts,counts)

      }

      if (level[i]=="genus")
      {
        joined_counts <- rbind(joined_counts,lv_g[,-1])
      }

      if (level[i]=="family")
      {
        joined_counts <- rbind(joined_counts,lv_f[,-1])
      }

      if (level[i]=="order")
      {
        joined_counts <- rbind(joined_counts,lv_o[,-1])
      }

      if (level[i]=="class")
      {
        joined_counts <- rbind(joined_counts,lv_c[,-1])
      }

      if (level[i]=="phylum")
      {
        joined_counts <- rbind(joined_counts,lv_p[,-1])
      }

    }


    joined_counts <- joined_counts[-1,]
    final_obj <- list(relabun=joined_counts,meta_data=covdata)
    return(final_obj)

  }


  ## If lowest level is genus

  if(length(chk[chk==5])==length(chk))

  {

    print("Lowest level is genus")

    # Clean OTU annotation and convert to composition abundance
    levels=as.data.frame(do.call('rbind',strsplit(rownames(counts),';')))

    if(missing(cleananno)){composition=FALSE}

    if(cleananno==TRUE)
    {
      levels$V6=gsub(pattern = '\\[',replacement = '',levels$V6)
      levels$V6=gsub(pattern = '\\]',replacement = '',levels$V6)
      levels=as.data.frame(apply(levels,2,function(x){
        y=gsub(' ','.',x)
        y=gsub('-','.',y)

        final_obj <- list(relabun=y,meta_data=covdata)
        return(final_obj)

      }
      ))}

    colnames(levels)=c('kingdom','phylum','class','order','family','genus')
    otu.names=paste('k__',levels[,1],'.p__',levels[,2],'.c__',levels[,3],'.o__',levels[,4],'.f__',levels[,5],'.g__',levels[,6],sep = '')
    rownames(counts)=otu.names

    if(composition==TRUE)
    {
      counts=t(t(counts)/colSums(counts))
    }



    ## Phylum

    lv_p=stats::aggregate(counts,by=list(levels$phylum),FUN=sum)
    order.id=order(levels$phylum)
    rownames(lv_p)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],sep = ''))

    ## Class
    lv_c=stats::aggregate(counts,by=list(levels$class),FUN=sum)
    order.id=order(levels$class)
    rownames(lv_c)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],sep = ''))

    ## Order
    lv_o=stats::aggregate(counts,by=list(levels$order),FUN=sum)
    order.id=order(levels$order)
    rownames(lv_o)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],sep = ''))

    ## Family
    lv_f=stats::aggregate(counts,by=list(levels$family),FUN=sum)
    order.id=order(levels$family)
    rownames(lv_f)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],'.f__',levels[order.id,5],sep = ''))
    ##Genus
    lv_g=stats::aggregate(counts,by=list(levels$genus),FUN=sum)
    order.id=order(levels$genus)
    rownames(lv_g)=unique(paste('k__',levels[order.id,1],'.p__',levels[order.id,2],'.c__',levels[order.id,3],
                                '.o__',levels[order.id,4],'.f__',levels[order.id,5],'.g__',levels[order.id,6],sep = ''))


    joined_counts <- matrix(nrow = 1,ncol = dim(lv_g)[2],dimnames = list(c(),c(colnames(lv_g))))

    for(i in 1:length(level))
    {

      # if (level[i]=="species")
      # {
      #   joined_counts <- rbind(joined_counts,counts)
      #
      # }

      if (level[i]=="genus")
      {
        joined_counts <- rbind(joined_counts,lv_g[,-1])
      }

      if (level[i]=="family")
      {
        joined_counts <- rbind(joined_counts,lv_f[,-1])
      }

      if (level[i]=="order")
      {
        joined_counts <- rbind(joined_counts,lv_o[,-1])
      }

      if (level[i]=="class")
      {
        joined_counts <- rbind(joined_counts,lv_c[,-1])
      }

      if (level[i]=="phylum")
      {
        joined_counts <- rbind(joined_counts,lv_p[,-1])
      }

    }

    joined_counts <- joined_counts[-1,]
    final_obj <- list(relabun=joined_counts,meta_data=covdata)
    return(final_obj)

  }

}





