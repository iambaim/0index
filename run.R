#install.packages(c("remotes", "RstoxData", "worms"))
#remotes::install_github("StoXProject/RstoxBase")

library(RstoxData)
library(openxlsx)
library(data.table)
library(worms)
library(sf)
library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)
library(scales)
library(gstat)
library(automap)
library(gridExtra)

defProj <- "+proj=longlat +datum=WGS84"

getDK <- function(fdm) {
    trTab <- c('9001'= 1,
            '9002' = 2,
            '9003' = 3,
            '9004' = 4,
            '9005' = 5,
            '9006' = 1,
            '9007' = 2,
            '9008' = 3,
            '9009' = 1,
            '9010' = 2,
            '9011' = 1,
            '9012' = 2,
            '9013' = 1,
            '9014' = 2,
            '9018' = 3,
            '9019' = 7
        )
    
    ret <- ifelse(fdm >= 9000, trTab[as.character(fdm)], floor(fdm/20))

    return(ret)
}

getScientificNames <- function(data_ref) {
    # Source Edvin's script
    source("https://github.com/Sea2Data/cruisetools/raw/master/taxaAnnotation/annotateTaxa.R")

    # Get taxa taxaTable
    ## Get list of aphias
    aphias <- unlist(unique(data_ref$catchsample[!is.na(aphia), "aphia"]))

    # Use cache
    cacheFile <- "taxaCache.rds"

    if(file.exists(cacheFile)) {
        print("Using taxa table cache...")
        taxaTable <- readRDS(cacheFile)
        # Just get the unidentified aphias
        aphias <- setdiff(aphias, taxaTable$AphiaID)
        if(length(aphias) > 0) {
            print(paste("Getting additional", length(aphias), "aphia IDs."))
            taxaTableNew <- makeTaxaTable(aphias)
            taxaTable <- rbind(taxaTable, taxaTableNew)
        }
    } else {
        taxaTable <- makeTaxaTable(aphias)
    }

    ## Save taxa table cache
    saveRDS(taxaTable, cacheFile)

    return(taxaTable)
}

processYear <- function(year, outdir = "results") {
    # Get the files
    files <- list.files(path = as.character(year), pattern="*.xml", recursive = TRUE, full.names = TRUE)

    print(files)

    # Out file
    dir.create(paste0("./", year, "/", outdir), recursive = TRUE, showWarnings = FALSE)
    outpath <- paste0("./", year, "/", outdir, "/", year, "_")

    # Read data
    raw <- ReadBiotic(files)

    # Table names
    wb_names <- c("mission", "fishstation", "catchsample", "individual")

    # Target columns
    targetColumns <- list(
        mission = c("startyear", "platform", "missionnumber", "platformname", "cruise"),
        fishstation = c("startyear", "platform", "missionnumber", "serialnumber", "nation", "catchplatform", "station",
            "fixedstation", "stationstartdate", "stationstarttime", "stationstopdate", "stationstoptime", "stationtype", "latitudestart",
            "longitudestart", "latitudeend", "longitudeend", "system", "area", "location", "bottomdepthstart", "bottomdepthstop", "bottomdepthmean",
            "fishingdepthstart", "fishingdepthstop", "fishingdepthcount", "fishingdepthmax", "fishingdepthmin", "fishingdepthmean",
            "fishingdepthtemperature", "gearno", "sweeplength", "gear", "gearcount", "direction", "gearflow", "vesselspeed", "logstart",
            "logstop", "distance", "gearcondition", "samplequality"),
        catchsample = c("startyear", "platform", "missionnumber", "serialnumber", "catchsampleid", "commonname", "catchcategory",
            "catchpartnumber", "aphia", "scientificname", "identification", "foreignobject", "sampletype", "group", "conservation",
            "catchproducttype", "raisingfactor", "catchweight", "catchvolume", "catchcount", "abundancecategory", "sampleproducttype",
            "lengthmeasurement", "lengthsampleweight", "lengthsamplevolume", "lengthsamplecount", "specimensamplecount", "agesamplecount"),
        individual = c("startyear", "platform", "missionnumber", "serialnumber", "catchsampleid", "specimenid", "individualproducttype",
            "individualweight", "individualvolume", "lengthresolution", "length")
    )

    # Deleted columns
    skipColumns <- lapply(wb_names, function(x) setdiff(names(raw[[1]][[x]]),targetColumns[[x]]))
    skipCol <- c()
    for(xx in seq_along(skipColumns)){
        skipCol <- union(skipCol, skipColumns[[xx]])
    }

    # Merge similar tables between files using data.table's rbindlist
    dt.u <- lapply(wb_names, function(x) rbindlist(lapply(raw, "[[", x)))
    names(dt.u) <- wb_names

    # Populate scientific name
    taxaTable <- getScientificNames(dt.u)
    dt.u$catchsample <- merge(dt.u$catchsample[,-c("scientificname")], taxaTable[, c("AphiaID", "scientificname")], 
                                by.x="aphia", by.y="AphiaID", all.x = TRUE)

    # Prepare Excel worbook
    wb <- createWorkbook()
    for(i in wb_names) {
        name <- paste0(i, "_unfiltered")
        addWorksheet(wb, name)
        writeData(wb, name, dt.u[[i]][, targetColumns[[i]], with = FALSE])
    }

    # Save workbook
    saveWorkbook(wb, file = paste0(outpath, "raw-unfiltered.xlsx"), overwrite = TRUE)


    # Prepare filter
    filt <- list(fishstation = c("stationtype %notin% c(2)",
                                    #"gear %in% c(3513)",
                                    "gear %in% c(3513, 3514)",
                                    "gearcondition %in% c(1)",
                                    "samplequality %in% c(1)"),
                catchsample = c("group %in% c(10, 13)")
    )
    filters <- lapply(raw, function(x) filt)

    # Filter the data
    biotic <- FilterBiotic(raw, filters)

    # Merge similar tables between files using data.table's rbindlist
    dt <- lapply(wb_names, function(x) rbindlist(lapply(biotic, "[[", x)))
    names(dt) <- wb_names
    dt$catchsample <- merge(dt$catchsample[,-c("scientificname")], taxaTable[, c("AphiaID", "scientificname")], 
                                by.x="aphia", by.y="AphiaID", all.x = TRUE)


    # Prepare Excel worbook
    wb <- createWorkbook()
    for(i in wb_names) {
        name <- i
        addWorksheet(wb, name)
        writeData(wb, name, dt[[i]][, targetColumns[[i]], with = FALSE])
    }

    # Save workbook
    saveWorkbook(wb, file = paste0(outpath, "raw-filtered.xlsx"), overwrite = TRUE)

    # merge all data/flattening all tables (all.x means include empty stations too)

    ## 1) Fishstations
    fishstations <- merge(dt[["mission"]], dt[["fishstation"]])

    ## Add stratum
    shp_xml <- read_xml("TIBIA_WGIBAR_polygoner.shp.xml")
    xml_proj <- read_xml(xml_text(xml_find_all(shp_xml, "//peXml")[[1]]))
    crs_wkt <- xml_text(xml_find_all(xml_proj, "//WKT")[[1]])
    StratumPolygon = st_read("TIBIA_WGIBAR_polygoner.shp")
    st_crs(StratumPolygon) <- crs_wkt

    # Remove stations without position
    #initialLen <- nrow(fishstations)
    #fishstations <- fishstations[!(is.na(latitudestart) | is.na(longitudestart)), ]
    #print(paste("Removing", initialLen - nrow(fishstations), "records from fishstation table due to missing coordinates!!!"))

    # Fill missing lat/lon
    fishstations[is.na(latitudeend), latitudeend:=latitudestart]
    fishstations[is.na(longitudeend), longitudeend:=longitudestart]


    fishstations[, lat:=(latitudestart + latitudeend) / 2]
    fishstations[, lon:=(longitudestart + longitudeend) / 2]

    pts <- fishstations[, c("platformname", "lon", "lat")]
    pts <- st_as_sf(pts, coords = c("lon", "lat"), crs = defProj)
    pts <- st_transform(pts, crs_wkt)

    fishstations[, stratum:=sapply(st_intersects(pts, StratumPolygon), function(z) if (length(z)==0) NA_integer_ else StratumPolygon$Area[z[1]])]

    ## Try to plot stratum and stations
    ## Get maxmin
    xmax <- pmax(st_bbox(StratumPolygon)["xmax"], st_bbox(pts)["xmax"])
    xmin <- pmin(st_bbox(StratumPolygon)["xmin"], st_bbox(pts)["xmin"])
    ymax <- pmax(st_bbox(StratumPolygon)["ymax"], st_bbox(pts)["ymax"])
    ymin <- pmin(st_bbox(StratumPolygon)["ymin"], st_bbox(pts)["ymin"])
    bmaxmin <- matrix(c(xmin, ymin, xmax, ymax), nrow = 2, ncol = 2)

    world <- ne_countries(scale = "medium", returnclass = "sf")

    stationmap <- ggplot(data = world) +
        geom_sf(size = 0.2, fill = "gray90") +
        geom_sf(data = StratumPolygon, alpha = 0, size = 0.2) + 
        geom_sf(data = pts, aes(color = platformname), size = 1.5, shape = 21, fill = "darkred") +
        coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs = st_crs(crs_wkt), expand = TRUE) + 
        theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25), 
            panel.background = element_rect(fill = "aliceblue")) + 
         ggtitle(paste(year, "- Station overview"))

    ggsave(paste0(outpath, "stationinfo.pdf"), stationmap)

    ## 2) Catchsamples
    catchsamples <- merge(fishstations, dt[["catchsample"]], by = intersect(names(dt[["fishstation"]]), names(dt[["catchsample"]])))

    ## Fix for some catchcount NAs but populated lengthsamplecount
    catchsamples[is.na(catchcount) & !is.na(lengthsamplecount), `:=`(catchcount = lengthsamplecount, catchweight = lengthsampleweight)]

    # Treat duplicate aphia with different catchpartnumber
    # First, get both the unique and duplicated rows
    dupGroup <- c(intersect(names(dt[["fishstation"]]), names(dt[["catchsample"]])), "aphia")
    duplicateRow <- which(duplicated(catchsamples, by = dupGroup))
    uniqueRow <- unique(catchsamples, by = dupGroup)[unique(catchsamples[duplicateRow, ..dupGroup], by = dupGroup), on = dupGroup]

    # Aggregate to sum the below columns from the duplicated rows
    sumCols <- c("catchweight", "catchvolume", "catchcount", "lengthsampleweight", "lengthsamplevolume", "lengthsamplecount", "specimensamplecount", "agesamplecount")
    cleanedSamples <- catchsamples[, lapply(.SD, sum), by=dupGroup, .SDcols=sumCols]
    cleanedSamples <- cleanedSamples[unique(catchsamples[duplicateRow, ..dupGroup], by = dupGroup), on = dupGroup]

    # Sanity check!
    if(nrow(uniqueRow) != nrow(cleanedSamples))
        stop("Number of duplicate rows and aggregate results seems off. Please double check!!!")

    # Helper function
    toCond <- function(row) {
        rown <- names(row)
        expr <- character(0)
        for(i in seq_along(row)) {
            expr <- c(expr, paste0(rown[i], " %in% c(", paste(row[[i]], collapse = ","), ")"))
        }
        return(paste(expr, collapse = " & "))
    }

    for(rw in seq_len(nrow(uniqueRow))) {
        row <- uniqueRow[rw]

        # Get the reference catchsampleid
        pickid <- row[["catchsampleid"]]

        # Fill the unique rows with the aggregated values
        row_1 <- row[, ..dupGroup]
        rowcond <- toCond(as.list(row_1))

        for(cl in sumCols) {
            catchsamples[eval(parse(text = paste(rowcond, "& catchsampleid %in%", pickid))), (cl):=cleanedSamples[eval(parse(text = rowcond))][[cl]]]
        }

        # Get all posible catchsampleid for the duplicate rows
        allid <- catchsamples[eval(parse(text = rowcond)),][["catchsampleid"]]

        # Construct individual sample adjustments
        otherid <- allid[!(allid %in% pickid)]
        conditions <- as.list(row[,..dupGroup])
        conditions[["catchsampleid"]] <- otherid
        conditions[["aphia"]] <- NULL
        conditions <- toCond(conditions)

        # Override all individual data with the corrected catchsample
        dt[["individual"]][eval(parse(text = conditions)), catchsampleid:=pickid]
    }

    # Remove duplicated catchsamples as we don't need it anymore
    catchsamples <- catchsamples[!duplicateRow,]

    # Fill up fishingdepthcount
    catchsamples[is.na(fishingdepthcount), fishingdepthcount := getDK(fishingdepthmax)]

    if(any(catchsamples$fishingdepthcount > 8)) {
        print("Re-check fishingdepthcount!!! Now we will just remove it...")
        print(catchsamples[fishingdepthcount > 8,])
        catchsamples <- catchsamples[!(fishingdepthcount > 8),]
    }

    # Calculate density
    catchsamples[, density:=((catchcount * 1852) / (20 * (distance / fishingdepthcount)))]

    # Calculate biomass
    catchsamples[, biomass:=((catchweight * 1852) / (20 * (distance / fishingdepthcount)))]

    # Check density
    if(any(is.na(catchsamples[["density"]]))) {
        print(catchsamples[is.na(density),])
        print("Found NAs in density. Must calculate it properly!")
    }

    # Casting density
    densityTable <- dcast(catchsamples, as.formula(paste0(paste(intersect(names(dt[["fishstation"]]), names(dt[["catchsample"]])), collapse = "+"), " ~ scientificname")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "density")

    # Casting biomass
    biomassTable <- dcast(catchsamples, as.formula(paste0(paste(intersect(names(dt[["fishstation"]]), names(dt[["catchsample"]])), collapse = "+"), " ~ scientificname")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "biomass")

    ## Save to PelagicData table
    PelagicData <- copy(catchsamples)
    suppressWarnings(PelagicData[, (skipCol):=NULL])
    densityTable <- merge(densityTable, fishstations, by=intersect(names(densityTable), names(fishstations)))
    suppressWarnings(densityTable[, (skipCol):=NULL])
    biomassTable <- merge(biomassTable, fishstations, by=intersect(names(biomassTable), names(fishstations)))
    suppressWarnings(biomassTable[, (skipCol):=NULL])
    wb <- createWorkbook()
    addWorksheet(wb, "PelagicData")
    writeData(wb, "PelagicData", PelagicData)
    addWorksheet(wb, "DensityPelagicData")
    writeData(wb, "DensityPelagicData", densityTable)
    addWorksheet(wb, "BiomassPelagicData")
    writeData(wb, "BiomassPelagicData", biomassTable)

    saveWorkbook(wb, file = paste0(outpath, "PelagicData.xlsx"), overwrite = TRUE)

    ## 3) Individuals

    # Merge
    individuals <- merge(catchsamples, dt[["individual"]], by = intersect(names(dt[["catchsample"]]), names(dt[["individual"]])))

    # Calculate lengthgroup (in cm)
    individuals[, lgr:=as.integer((length*100)/0.5)*0.5+0.5/2]

    # Construct lengthgroup string
    getLGstring <- function(lgr) {
        base <- (lgr-0.25) * 10
        return(ifelse(is.na(base), NA, paste0(base, " - ", base + 4, " mm")))
    }
    individuals[, lengthgroup:=getLGstring(lgr)]

    ## Calculate keff
    individuals[, keff:=1]
    individuals[aphia %in% c(126436, 126437, 126441, 126433) & lgr < 14.7, keff:=17.065 * exp(-0.1932 * lgr)]
    individuals[aphia %in% c(126436, 126437, 126441, 126433) & lgr < 2.8, keff:=10]
    individuals[aphia %in% c(126417) & lgr < 9.8, keff:=357.23 * exp(-0.6007 * lgr)]
    individuals[aphia %in% c(126417) & lgr < 4.1, keff:=30]
    individuals[aphia %in% c(126735) & lgr < 11.7, keff:=7.2075 * exp(-0.1688 * lgr)]

    # Calculate fishlength per-station
    stratumGroup <- c("stratum", "aphia", "scientificname")
    stratumIndLength <- individuals[, list(avg = mean(length), min = min(length), max = max(length), keff = sum(keff)), by = stratumGroup]

    # Calculate fishlength in the whole sea
    wholeGroup <- c("aphia", "scientificname")
    wholeIndLength <- individuals[, list(avg = mean(length), min = min(length), max = max(length), keff = sum(keff)), by = wholeGroup]

    # Calculate fishlength
    individualGroup <- c("platformname", "station", "stratum", "serialnumber", "stationstartdate", "stationstarttime", "stationstopdate", "stationstoptime", "lat", "lon", "distance", "fishingdepthcount", "aphia", "scientificname", "catchcount", "lengthsamplecount", "lgr", "lengthgroup")
    indSummary <- individuals[, list(avg = mean(length), min = min(length), max = max(length), N = .N, proportion = .N/lengthsamplecount, keff = sum(keff)), by = individualGroup]

    ## Calculate length distribution
    indSummary[, lengthdistribution:=catchcount * proportion]
    indSummary[, lengthdistributionKeff:=catchcount/N * keff * proportion]

    # Prepare wide summary for per-length-group
    individualGroup_1 <- individualGroup[!individualGroup %in% c("lgr", "lengthgroup")]
    widesum1 <- dcast(indSummary, as.formula(paste0(paste(individualGroup_1, collapse = "+"), " ~ lengthgroup")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "N")
    widesum2 <- dcast(indSummary, as.formula(paste0(paste(individualGroup_1, collapse = "+"), " ~ lengthgroup")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "proportion")
    widesum3 <- dcast(indSummary, as.formula(paste0(paste(individualGroup_1, collapse = "+"), " ~ lengthgroup")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "lengthdistribution")
    widesum4 <- dcast(indSummary, as.formula(paste0(paste(individualGroup_1, collapse = "+"), " ~ lengthgroup")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "keff")
    widesum5 <- dcast(indSummary, as.formula(paste0(paste(individualGroup_1, collapse = "+"), " ~ lengthgroup")), fun.agg = function(x)  sum(x, na.rm = TRUE), value.var = "lengthdistributionKeff")


    # Post-processing keff summary
    sumKeff <- copy(widesum5)
    mmCols <- colnames(sumKeff)[grep("mm+",colnames(sumKeff))]
    sumKeff[, sumkeff := rowSums(.SD), .SDcols = mmCols]
    sumKeff <- sumKeff[aphia %in% c(126436, 126437, 126441, 126433, 126417, 126735),]
    sumKeff[, density:=((catchcount * 1852) / (20 * (distance / fishingdepthcount)))]
    sumKeff[, densityKeff:=((sumkeff * 1852) / (20 * (distance / fishingdepthcount)))]

    # Combined data
    suppressWarnings(indSummary[, (skipCol):=NULL])
    suppressWarnings(widesum1[, (skipCol):=NULL])
    suppressWarnings(widesum2[, (skipCol):=NULL])
    suppressWarnings(widesum3[, (skipCol):=NULL])
    suppressWarnings(widesum4[, (skipCol):=NULL])
    suppressWarnings(widesum5[, (skipCol):=NULL])
    sumKeff_1 <- sumKeff
    suppressWarnings(sumKeff_1[, (c(skipCol, mmCols)):=NULL])

    # Prettify column ordering
    lengths <- colnames(widesum1)[grep("mm", colnames(widesum1))]
    lengths <- lengths[order(as.numeric(lapply(strsplit(lengths, " "), "[[", 1)))]
    newcols <- c(colnames(widesum1)[-grep("mm", colnames(widesum1))], lengths)

    # Histograms
    LengthFreq <- widesum3[, lapply(.SD, sum), by=scientificname, .SDcols=(lengths)]
    LengthFreqKeff <- widesum5[aphia %in% c(126436, 126437, 126441, 126433, 126417, 126735), lapply(.SD, sum), by=scientificname, .SDcols=(lengths)]

    histo <- LengthFreq[, c(list(scientificname = scientificname), lapply(.SD, function(x) (x/rowSums(.SD))*100)), .SDcols = lengths]
    histoKeff <- LengthFreqKeff[, c(list(scientificname = scientificname), lapply(.SD, function(x) (x/rowSums(.SD))*100)), .SDcols = lengths]

    setcolorder(widesum1, newcols)
    setcolorder(widesum2, newcols)
    setcolorder(widesum3, newcols)
    setcolorder(widesum4, newcols)
    setcolorder(widesum5, newcols)

    # Make summary
    wb <- createWorkbook()

    addWorksheet(wb, "StratumIndLength")
    addWorksheet(wb, "WholeIndLength")

    addWorksheet(wb, "summary")
    addWorksheet(wb, "summaryWide1-N")
    addWorksheet(wb, "summaryWide2-prop")
    addWorksheet(wb, "summaryWide3-lengthdist")
    addWorksheet(wb, "summaryWide4-keff")
    addWorksheet(wb, "summaryWide5-lengthdistkeff")
    addWorksheet(wb, "summary-keff")

    writeData(wb, "StratumIndLength", stratumIndLength)
    writeData(wb, "WholeIndLength", wholeIndLength)

    writeData(wb, "summary", indSummary)
    writeData(wb, "summaryWide1-N", widesum1)
    writeData(wb, "summaryWide2-prop", widesum2)
    writeData(wb, "summaryWide3-lengthdist", widesum3)
    writeData(wb, "summaryWide4-keff", widesum4)
    writeData(wb, "summaryWide5-lengthdistkeff", widesum5)
    writeData(wb, "summary-keff", sumKeff_1)

    saveWorkbook(wb, file = paste0(outpath, "summary.xlsx"), overwrite = TRUE)

    # Create histograms
    wb <- createWorkbook()

    addWorksheet(wb, "Length frequency")
    addWorksheet(wb, "Length frequency Keff")
    writeData(wb, "Length frequency", t(histo), colNames = FALSE, rowNames = TRUE)
    writeData(wb, "Length frequency Keff", t(histoKeff), colNames = FALSE, rowNames = TRUE)

    saveWorkbook(wb, file = paste0(outpath, "LengthFrequency.xlsx"), overwrite = TRUE)

    # Prepare grid
    krig.grid <- st_sample(StratumPolygon, 10000, type="regular")

    plotKrig <- function (ap, krig.grid, sumKeff, shp, yr) {

        mybreaks <- 10^c(2:6)

        keffmap <- sumKeff[scientificname == ap,]
        keffmap <- st_as_sf(keffmap, coords = c("lon", "lat"), crs = defProj)
        keffmap <- st_transform(keffmap, crs_wkt)

        plot.density <- ggplot(data = world) +
            geom_sf(size = 0.2, fill = "gray90") +
            geom_sf(data = keffmap, aes(color = densityKeff, size = densityKeff), alpha=0.5) + 
            geom_sf(data = shp, fill='transparent', color = "black", size = 0.2) +
            coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs=st_crs(crs_wkt), expand = TRUE) +
            scale_size("Density", trans="log10", breaks=mybreaks) +
            scale_color_binned("Density", trans="log10", breaks=mybreaks, type = "viridis") +
            guides(color=guide_legend(), size = guide_legend()) +
            theme_bw() + ggtitle(paste(yr, ap, "- Density w/ Keff"))

        print(ap)
        print(nrow(keffmap))

        kriging_result <- tryCatch(
                                {
                                    autoKrige(log10(densityKeff)~1, as(keffmap, "Spatial"), as(krig.grid, "Spatial"))
                                },
                                error=function(cond) {
                                    message(cond)
                                    return(NA)
                                })
        if(!is.na(kriging_result[[1]])) {
            dt.krig <- as.data.table(kriging_result$krige_output)
            plot.density.pred <- ggplot(data = world) +
                geom_sf(size = 0.2, fill = "gray90") +
                geom_raster(data=dt.krig, aes(x=coords.x1, y=coords.x2, fill=var1.pred))+
                geom_sf(data = shp, fill='transparent', color = "black", size = 0.2) +
                coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs = st_crs(crs_wkt), expand = TRUE) +
                scale_fill_gradientn("Density\n(log10)", breaks=log10(mybreaks), colours = c("orange", "yellow", "green",  "sky blue", "blue")) + 
                theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ggtitle(paste(yr, ap, "- Density w/ Keff (prediction)"))
        } else {
            plot.density.pred <- ggplot(data = world) +
                geom_sf(size = 0.2, fill = "gray90") +
                geom_sf(data = shp, fill='transparent', color = "black", size = 0.2) +
                coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs = st_crs(crs_wkt), expand = TRUE) +
                theme_bw() + ggtitle(paste(yr, ap, "w/ Keff (Kriging error!!!)"))
        }

        return(list(plot.density, plot.density.pred))
    }

    species <- unique(sumKeff$scientificname)
    plots <- lapply(species, plotKrig, krig.grid, sumKeff, StratumPolygon, year)

    ggsave(
        filename = paste0(outpath, "density-keff-plots.pdf"), 
        plot = marrangeGrob(unlist(plots, recursive = FALSE), nrow=1, ncol=2, top = ""),
        width = 15, height = 9
    )

    mybreaks <- 10^c(2:6)

    keffData <- indSummary[aphia %in% c(126436, 126437, 126441, 126433, 126417, 126735),]

    kData <- keffData[, c("lon", "lat", "scientificname", "lengthdistributionKeff", "N", "lgr")]
    kData <- st_as_sf(kData, coords = c("lon", "lat"), crs = defProj)
    kData <- st_transform(kData, crs_wkt)

    plots <- list()

    for(sp in unique(kData$scientificname)) {

        keffPlot <- ggplot(data = world) +
            geom_sf(size = 0.2, fill = "gray90") + 
            #facet_wrap(vars(scientificname), nrow = 2) +
            geom_sf(data = kData[kData$scientificname == sp, ], aes(size=lengthdistributionKeff, color=lgr), alpha=0.5) +
            geom_sf(data = StratumPolygon, fill='transparent', color = "black", size = 0.2) +
            coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs = st_crs(crs_wkt), expand = TRUE) +
            scale_size("Length Distribution") +
            scale_colour_gradientn("Length Group", colours = hcl.colors(13, "Dark 3")) +
            theme_bw() + theme(legend.position="top") + ggtitle(paste(year, sp, "- Length Distribution"))

        lengthPlot <- ggplot(data = world) +
            geom_sf(size = 0.2, fill = "gray90") +
            #facet_wrap(vars(scientificname), nrow = 2) +
            geom_sf(data = kData[kData$scientificname == sp, ], aes(size=N, color=lgr), alpha=0.5) +
            geom_sf(data = StratumPolygon, fill='transparent', color = "black", size = 0.2) +
            coord_sf(xlim = bmaxmin[1,], ylim = bmaxmin[2,], crs = st_crs(crs_wkt), expand = TRUE) +
            scale_size("N") +
            scale_colour_gradientn("Length Group", breaks = c(1:15), colours = hcl.colors(15, "Dark 3")) +
            theme_bw() + theme(legend.position="top") + ggtitle(paste(year, sp, "- Number of Sample by Length Groups"))

        plots <- c(plots, list(keffPlot, lengthPlot))

    }

    ggsave(
        filename = paste0(outpath, "lengthdist-plots.pdf"), 
        plot = marrangeGrob(plots, nrow=1, ncol=2, top = ""),
        width = 15, height = 9
    )
}



