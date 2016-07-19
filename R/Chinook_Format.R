#' @title Format Chinook Salmon JTRAP data needed for MYTSBE models
#'
#' @description This function reconfigures Chinook Salmon JTRAP data into the appropriate capture-mark-recapture format needed for the Multi-Year Time Stratified Bayesian Estimator (MYTSBE) models.
#'
#' Chinook "TU" dispositioned fish are marked individuals released upstream (nraw) of the rotary screw trap (RST) on a given date. "RE RC" dispositioned Chinook are recaptured individuals (m) at the RST on a given date.
#' The sum of "TU", "TD", "NTT", "NTD", "NTR" & "NTS" dispositioned fish are captured unmarked individuals (u) at the RST on a given date.
#'
#' Fry and NTS disposition individuals captured between the earliest trapping date to the specified smolt cut-off are broken into "YOYnraw", "YOYm", and "YOYu". Precocial Chinook, designated with precocial disposition, are excluded from the formatting. Yearling fish migrating after the specified smolt date are excluded from the formatting.
#'
#' Years are standardized in that each year begins at the earliest trapping date and ends at the latest trapping date the RST was ever in operation since RST installation.
#'
#'  PTAGIS data are needed to calculate the likelihood an individual is available for recapture during a proceeding day.
#'  Movement parameters of individauls are calculated using a MLE approach beginning with Day 1 (day after release) through Day 5, for the entirety of the data set. After day 5, it is assumed that migrating fish have traveled passed the RST
#' or left the population. "n" is the estimated number of fish available for recapture on a given date rounded to the nearest integer.
#'
#' This package exports a .csv file to the working directory titled "Trap name, species, today's date, formatted data.csv"
#'
#' @param JTRAP_data JTRAP data query
#' @param PTAG_data PTAGIS data query (needs to be specific to the species of interest)
#' @param strata length of desired strata in days (e.g. 1 week -> 7 days)
#' @param species character string needed to subset data (should be identical to the species name in the JTRAP_data set)
#' @param trap.name character string used for titles and descriptions of reports
#' @param smolt.date "MM-DD" smolt classification date needed for formatting
#' @import plyr
#' @importFrom useful shift.column
#' @export
#' @return NULL

Format_Chinook <- function(JTRAP_data,
                         PTAG_data,
                         strata = 7,
                         smolt.date = "07-01",
                         species,
                         trap.name){

  JTRAP_data$PT2Date <-as.Date(JTRAP_data$PT2Date, format = "%m/%d/%Y")
  JTRAP_data$Year <-as.numeric(format(JTRAP_data$PT2Date, format = "%Y"))
  JTRAP_data$jday <- yday(JTRAP_data$PT2Date)
  JTRAP_data$Smolt_pot <- yday(as.Date(paste(JTRAP_data$Year,"-",smolt.date, sep ="")))
  JTRAP_data$BroodYr <- as.numeric(JTRAP_data$BroodYr)
  JTRAP_data$NFish <- as.numeric(JTRAP_data$NFish)


  day.effort <- count(JTRAP_data, "PT2Date")
  day.effort$freq <- 1
  colnames(day.effort)<- c("date", "effort")


  # movement parameter from PTAGIS query
  PTAG_data$Days.Between <- as.numeric(as.Date(PTAG_data$Recap.Date.MMDDYYYY, format = "%m/%d/%Y") - as.Date(PTAG_data$Mark.Date.MMDDYYYY,format = "%m/%d/%Y")) #find days between mark and recapture
  PTAG_data <- PTAG_data[PTAG_data$Days.Between > 0 ,] # clean data with negative days between mark/recapture
  PTAG_data <- PTAG_data[PTAG_data$Days.Between <= 10 ,] # clean data by removing recapture days over 10
  PTAG_data <- PTAG_data[order(PTAG_data$Tag.Code, PTAG_data$Days.Between, decreasing = FALSE),] # order by Tag.Code & Days.Between so the shortest "Days.Between" values are retained for duplicates. These fish are presumed to be recaptured fish released upstream and caught again recaptured
  PTAG_data <- PTAG_data[!duplicated(PTAG_data$Tag.Code),] # remove duplicate PTAG codes

  Mig_days <- count(PTAG_data, "Days.Between") #count total # of fish

  Mig_days$Percent <- Mig_days$freq/sum(Mig_days$freq) #Create percentages

  ################################
  #        nraw,m,u formatting      #
  ################################
  # Choose  species dataset from original raw input
  JTRAP_data <- JTRAP_data[JTRAP_data$SpName == species , ]

  JTRAP_fry <- subset(JTRAP_data, ((JTRAP_data$jday < JTRAP_data$Smolt_pot) &
                                     ((JTRAP_data$Disposition_Acronym == "NTS") |
                                        (JTRAP_data$Disposition_Acronym == "BBY") |
                                        (JTRAP_data$Disposition_Acronym == "RE BY") |
                                        (JTRAP_data$Disposition_Acronym == "YOY")
                                     )))

  JTRAP_data <- subset(JTRAP_data, !((JTRAP_data$jday < JTRAP_data$Smolt_pot) &
                                       ((JTRAP_data$Disposition_Acronym == "NTS") |
                                          (JTRAP_data$Disposition_Acronym == "BBY") |
                                          (JTRAP_data$Disposition_Acronym == "RE BY") |
                                          (JTRAP_data$Disposition_Acronym == "YOY") |
                                          ((JTRAP_data$Year-JTRAP_data$BroodYr) == 1))))

  trap_m <- JTRAP_data[JTRAP_data$Disposition_Acronym =="RE RC" ,]
  trap_n <- JTRAP_data[JTRAP_data$Disposition_Acronym =="TU" ,]
  trap_u <- JTRAP_data[(JTRAP_data$Disposition_Acronym == "TU") | (JTRAP_data$Disposition_Acronym == "TD") | (JTRAP_data$Disposition_Acronym == "NTT") | (JTRAP_data$Disposition_Acronym == "NTD") | (JTRAP_data$Disposition_Acronym == "NTR") | (JTRAP_data$Disposition_Acronym == "NTS"), ]

  # Count the number of times the "PT2Date" shows up and weight it by "NFish" to get daily counts
  trap_m_total <- count(trap_m, "PT2Date", "NFish")
  trap_n_total <- count(trap_n, "PT2Date", "NFish")
  trap_u_total <- count(trap_u, "PT2Date", "NFish")

  # Merge the total count data sets back together
  trap_merge <- merge(trap_m_total, trap_n_total, by="PT2Date", all=TRUE)
  colnames(trap_merge)<- c("PT2Date","m","nraw")
  trap_merge <- merge(trap_merge, trap_u_total, by="PT2Date" , all=TRUE)
  colnames(trap_merge)[4]<- "u"

  names(trap_merge) <- c("date","m","nraw","u")
  trap_merge$date <- as.Date(trap_merge$date, format = "%m/%d/%Y")

  minyear <- as.numeric(min(format(trap_merge$date, format= "%Y"), na.rm=TRUE))
  maxyear <- as.numeric(max(format(trap_merge$date, format= "%Y"), na.rm=TRUE))

  dates <- as.data.frame(seq(from = as.Date(paste(minyear,"-01-01",sep="")),
                             to = as.Date(paste(maxyear,"-12-31",sep="")),
                             by = "days"))
  names(dates) <- c("date")

  ###############################################
  #        Calendar fill, days, and effort      #
  ###############################################

  # merge with calendar with trap data
  trap_date_filled = merge(trap_merge, dates, by="date", all=TRUE)

  ###############################################
  #    Adding Julian Date & Year Variables      #
  ###############################################

  trap_date_filled$julian = yday(trap_date_filled$date)
  trap_date_filled$year = format(trap_date_filled$date, "%Y")

  trap_full <-trap_date_filled

  #########################################
  #         Incorporating travel time     #
  #########################################

  # create matrix when a fish released at time 0 will be available for recapture
  trap_mig <- as.data.frame(trap_full$nraw*Mig_days[1,3])
  trap_mig$Nd2<- trap_full$nraw*Mig_days[2,3]
  trap_mig$Nd3<- trap_full$nraw*Mig_days[3,3]
  trap_mig$Nd4<- trap_full$nraw*Mig_days[4,3]
  trap_mig$Nd5<- trap_full$nraw*Mig_days[5,3]

  #apply a function that sums the diagonal of the matrix from lower left to upper right for each day
  trap_mig <- as.matrix(trap_mig)
  trap_mig[is.na(trap_mig)] <- 0
  vals <- as.data.frame(sapply(split(as.matrix(trap_mig), row(trap_mig) + col(trap_mig)), sum))
  #round the variables
  vals <- round(vals)
  #rename the header
  colnames(vals) <- c("n")
  #remove excess rows
  vals <- head(vals, -4)
  #add adjusted n back to trap data
  trap_full$n <- vals$n

  #shift n down one to the appropriate day
  trap_full <- shift.column(trap_full, "n", len=1, up=FALSE)
  trap_full$n <- trap_full$n.Shifted
  trap_full$n.Shifted <- NULL

  trap_full$n <- ifelse((trap_full$m > 0) | (trap_full$nraw > 0) | (trap_full$u > 0), trap_full$n, 0)

  #########################################
  #                YOY                    #
  #########################################

  trap_yoym <- JTRAP_fry[JTRAP_fry$Disposition_Acronym =="RE BY" ,]
  trap_yoyn <- JTRAP_fry[JTRAP_fry$Disposition_Acronym =="BBY" ,]
  trap_yoyu <- JTRAP_fry[(JTRAP_fry$Disposition_Acronym == "YOY") | (JTRAP_fry$Disposition_Acronym == "BBY") | (JTRAP_fry$Disposition_Acronym == "NTS"), ]

  # Count the number of times the "PT2Date" shows up and weight it by "NFish" to get daily counts
  trap_yoym_total <- count(trap_yoym, "PT2Date", "NFish")
  trap_yoyn_total <- count(trap_yoyn, "PT2Date", "NFish")
  trap_yoyu_total <- count(trap_yoyu, "PT2Date", "NFish")

  # Merge the total count data sets back together
  trap_merge_yoy <- merge(trap_yoym_total, trap_yoyn_total, by="PT2Date", all=TRUE)
  colnames(trap_merge)<- c("PT2Date","yoym","yoynraw")
  trap_merge_yoy <- merge(trap_merge_yoy, trap_yoyu_total, by="PT2Date" , all=TRUE)
  colnames(trap_merge_yoy)[4]<- "yoyu"

  names(trap_merge_yoy) <- c("date","yoym","yoynraw","yoyu")
  trap_merge_yoy$date <- as.Date(trap_merge_yoy$date, format = "%m/%d/%Y")

  trap_date_yoy_filled = merge(trap_merge_yoy, dates, by="date", all=TRUE)

  trap_date_yoy_filled <- shift.column(trap_date_yoy_filled, "yoynraw", len=1, up=FALSE)
  trap_date_yoy_filled$yoyn <- trap_date_yoy_filled$yoynraw.Shifted
  trap_date_yoy_filled$yoynraw.Shifted <- NULL

  trap_full <- merge(trap_full, trap_date_yoy_filled, by="date", all=TRUE)

  #add a column "days" and fill with "1"s
  trap_full["days"] <- 1

  #create an "effort" column that is filled with 1 if a fish was captured or a 0 if a fish wasn't captured
  trap_full <- merge(trap_full, day.effort, by="date", all=TRUE)

  #if the trap effort is "NA", change it to "0"
  trap_full$effort[is.na (trap_full$effort)] = 0

  #########################################
  #              Stratafying              #
  #########################################

  trap_full$strata=cut(trap_full$julian, seq(0, max(trap_full$julian, na.rm = TRUE)+7, by=strata), labels=FALSE)

  #Summarize by Year and strata for Chinook
  trap_week=ddply(trap_full, c("year","strata"), summarize,
                  m=sum(m, na.rm = TRUE),
                  nraw=sum(nraw, na.rm = TRUE),
                  u=sum(u, na.rm = TRUE),
                  n=sum(n,na.rm=TRUE),
                  yoym=sum(yoym, na.rm = TRUE),
                  yoynraw=sum(yoynraw, na.rm = TRUE),
                  yoyu=sum(yoyu, na.rm = TRUE),
                  yoyn=sum(yoyn, na.rm = TRUE),
                  days=sum(days, na.rm = TRUE),
                  effort=sum(effort, na.rm = TRUE)
  )

  #Find the first and last Jweek the trap has operated through
  #jweek
  weeks=count(trap_week, "strata", "effort")
  trap_season=weeks[weeks$freq > 0 ,]
  trap_start=as.numeric(min(trap_season$strata, na.rm=TRUE))
  trap_finish=as.numeric(max(trap_season$strata, na.rm=TRUE))

  #Find the first year the trap was in operation
  years=count(trap_week, "year", "effort")
  trap_year=years[years$freq > 0 ,]
  trap_start_year=as.numeric(min(trap_year$year, na.rm=TRUE))

  #Chop trap data to only the 1st year of operation, first ordinal week, and last ordinal week
  trap_complete=trap_week[trap_week$strata >= trap_start , ]
  trap_complete=trap_complete[trap_complete$strata <= trap_finish ,]
  trap_complete=trap_complete[trap_complete$year >= trap_start_year ,]

  trap_complete$u<- ifelse((trap_complete$effort > 0) , trap_complete$u, NA)
  trap_complete$m<- ifelse((trap_complete$m > trap_complete$n), trap_complete$n, trap_complete$m)
  trap_complete$yoyu<- ifelse((trap_complete$effort > 0) , trap_complete$yoyu, NA)
  trap_complete$yoym<- ifelse((trap_complete$yoym > trap_complete$yoyn), trap_complete$yoyn, trap_complete$yoym)
  trap_complete<- trap_complete[-nrow(trap_complete),]

  today <- Sys.Date()

  write.csv(trap_complete, file = paste(trap.name, species, today, "formatted data.csv"))
}
