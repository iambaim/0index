# Install requirements
requiredPKGS <- c("xml2", "RstoxData", "openxlsx", "data.table", "worms", "sf",
    "ggplot2", "rnaturalearthdata", "rnaturalearth", "scales",
    "gstat", "automap", "gridExtra")
newPKGS <- requiredPKGS[!requiredPKGS %in% installed.packages()[,"Package"]]
installStatus <- lapply(newPKGS, install.packages, repos="https://cloud.r-project.org/")

source("run.R")

# Full version: https://github.com/REDUS-IMR/stsdownloader/blob/master/R/utils.R
getBioticFiles <- function(cruiseNo, shipName, snapshot = NA, targetDir = "./") {
 
    # Fix ship name
    shipName <- gsub("\\s+\\([0-9]+\\)", "", shipName)
    
    # Fix for GO SARS
    if (shipName == "G.O. Sars")
        shipName <- "G.O.Sars"

    url <- paste0("http://tomcat7.imr.no:8080/apis/nmdapi/biotic/v3?type=findByCruise&shipname=", shipName, "&cruisenr=", cruiseNo)
    print(url)
    doc <- tryCatch({
        read_xml(URLencode(url))
    }, error = function(e) {
        NULL
    })

    if(is.null(doc)) {
        return(list(filename=paste0("biotic_cruiseNumber_", cruiseNo, "_", shipName, ".xml"), status = FALSE))
    }

    path <- xml_text(xml_find_all(doc, "//*[local-name() = 'element'][@name='path']"))
    print(path)

    if (is.na(snapshot)) {
        url <- paste0("http://tomcat7.imr.no:8080/apis/nmdapi/biotic/v3/", path, "/snapshot?version=3.0")
        print(url)
        doc <- read_xml(URLencode(url))

        snapshots <- xml_text(xml_find_all(doc, "//*[local-name() = 'element'][@name='snapshot time']"))
        print(snapshots)
        snapshots <- snapshots[snapshots != "latest"]
        snapshots <- snapshots[order(as.Date(snapshots), decreasing = TRUE)]
        print(snapshots)
        snapshot <- snapshots[1]
    }

     # Construct the correct file name with snapshot
    correctFile <- paste0("biotic_cruiseNumber_", cruiseNo, "_", shipName, "_", snapshot, ".xml")

    # Download latest file
    url <- paste0("http://tomcat7.imr.no:8080/apis/nmdapi/biotic/v3/", path, "/snapshot/", snapshot, "?version=3.0")

    fullOutDir <- paste0(targetDir, "/input/biotic")

    if(!dir.exists(fullOutDir))
        dir.create(fullOutDir, recursive = TRUE)

    status <- tryCatch({
        download.file(URLencode(url), paste0(fullOutDir, "/", correctFile))
    }, error = function(e) {
        FALSE
    })

    if(status == 0)
        status <- TRUE
    else
        status <- FALSE

    return(list(filename=correctFile, status = status))
}

fileTbl <- as.data.table(read.xlsx("filelist.xlsx"))

allRes <- data.table()

for(yr in "2020") {
    col <- colnames(fileTbl)
    subYr <- fileTbl[fileTbl[[col[1]]] == yr,]
    results <- apply(subYr, 1, function(x) getBioticFiles(x[[2]], x[[3]], targetDir = x[[1]]))
    tblRes <- rbindlist(results)
    tblRes[, year:=yr]
    allRes <- rbind(allRes, tblRes)

    processYear(yr)
}

# Save summary
wb <- createWorkbook()
addWorksheet(wb, "status")
writeData(wb, "status", allRes)
saveWorkbook(wb, file = "status.xlsx", overwrite = TRUE)
