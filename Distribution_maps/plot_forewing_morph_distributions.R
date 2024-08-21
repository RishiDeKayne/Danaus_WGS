

library(ggplot2)
library(maps)
library(mapproj)
library(geosphere)
library(scatterpie)

#a function that divides the world into blocks of a given size
#it returns the longitude and latitude 'name' of each block and its midpoint coordinates
define_all_blocks <- function(blocksize){
    #use the cut() and levels() functions to divide the world and get the 'name' of each division
    longitude_level_names <- levels(cut(0, breaks = seq(-180, 180, blocksize)))
    latitude_level_names <- levels(cut(0, breaks = seq(-90, 90, blocksize)))
    
    #get the longitude and latitude midpoint of each level
    midpoint_long <- seq(-180+blocksize/2, 180, blocksize)
    midpoint_lat <- seq(-90+blocksize/2, 90, blocksize)
    
    #and add the level "name" for each midpoint value
    names(midpoint_long) <- levels(longitude_level_names)
    names(midpoint_lat) <- levels(latitude_level_names)
    
    #now define ALL bloicks by taking all possible combinations of longitude and latitude
    blocks <- expand.grid(longitude_level_names, latitude_level_names)
    
    #add column names for the blocks
    names(blocks) <- c("block_name_long", "block_name_lat")
    
    #now add the midpoints for each block
    #this is added by retrieving the midpoint value using the block's longitude and latitude names    
    blocks$block_midpoint_long <- midpoint_long[blocks$block_name_long]
    blocks$block_midpoint_lat <- midpoint_lat[blocks$block_name_lat]
    
    blocks
    }

# a function that assigns each record to a block based on its coordinates
# it returns a list of rows giving the records that lie within each block
# many of the blocks will have zero records, because this outputs ALL blocks in the world,
# and many will be in the ocean etc.
group_records_into_blocks <- function(longitude, latitude, blocksize, block_data){
    #use cut() to assign a longitude and latitude level for each record
    block_assignment_long <- cut(longitude, breaks = seq(-180, 180, blocksize))
    block_assignment_lat <- cut(latitude, breaks = seq(-90, 90, blocksize))
    
    #now for each block in the world, retrieve the records that match the correct name for BOTH longitude AND latitude
    lapply(1:nrow(block_data), function(i) which(block_assignment_long == block_data[i,"block_name_long"] &
                                                 block_assignment_lat == block_data[i,"block_name_lat"]))
    }


cap <- function(x, max) ifelse(x <= max, x, max)

jitter <- function(x, dist) sapply(x, function(y) runif(1, y-dist, y+dist))




################# basic map #####################


# regions to include on map
regions = c("Afghanistan"  , "Angola"  , "Albania"  , "Finland"  , "Andorra"  , "United Arab Emirates"  , "Armenia"  , "American Samoa"  , "Australia"  , "French Southern and Antarctic Lands" , "Austria"  , "Azerbaijan"  , "Burundi"  , "Belgium"  , "Benin"  , "Burkina Faso"  , "Bangladesh"  , "Bulgaria"  , "Bahrain"  , "Bosnia and Herzegovina"  , "Belarus"  , "Brunei"  , "Bhutan"  , "Botswana"  , "Central African Republic"  , "Switzerland"  , "China"  , "Ivory Coast"  , "Cameroon"  , "Democratic Republic of the Congo"  , "Republic of Congo"  , "Cook Islands"  , "Comoros"  , "Cape Verde"  , "Cyprus"  , "Czech Republic"  , "Germany"  , "Djibouti"  , "Denmark"  , "Algeria"  , "Egypt"  , "Eritrea"  , "Canary Islands"  , "Spain"  , "Estonia"  , "Ethiopia"  , "Fiji"  , "Reunion"  , "Mayotte"  , "France"  , "Faroe Islands"  , "Micronesia"  , "Gabon"  , "UK"  , "Georgia"  , "Guernsey"  , "Ghana"  , "Guinea"  , "Gambia"  , "Guinea-Bissau"  , "Equatorial Guinea"  , "Greece"  , "Heard Island"  , "Croatia"  , "Hungary"  , "Indonesia"  , "Isle of Man"  , "India"  , "Cocos Islands"  , "Christmas Island"  , "Chagos Archipelago"  , "Ireland"  , "Iran"  , "Iraq"  , "Iceland"  , "Israel"  , "Italy"  , "San Marino"  , "Jersey"  , "Jordan"  , "Japan"  , "Siachen Glacier"  , "Kazakhstan"  , "Kenya"  , "Kyrgyzstan"  , "Cambodia"  , "Kiribati"  , "South Korea"  , "Kosovo"  , "Kuwait"  , "Laos"  , "Lebanon"  , "Liberia"  , "Libya"  , "Liechtenstein"  , "Sri Lanka"  , "Lesotho"  , "Lithuania"  , "Luxembourg"  , "Latvia"  , "Morocco"  , "Monaco"  , "Moldova"  , "Madagascar"  , "Maldives"  , "Marshall Islands"  , "Macedonia"  , "Mali"  , "Malta"  , "Myanmar"  , "Montenegro"  , "Mongolia"  , "Northern Mariana Islands"  , "Mozambique"  , "Mauritania"  , "Mauritius"  , "Malawi"  , "Malaysia"  , "Namibia"  , "New Caledonia"  , "Niger"  , "Norfolk Island"  , "Nigeria"  , "Niue"  , "Bonaire"  , "Netherlands"  , "Norway"  , "Nepal"  , "Nauru"  , "New Zealand"  , "Oman"  , "Pakistan"  , "Panama"  , "Pitcairn Islands"  , "Philippines"  , "Palau"  , "Papua New Guinea"  , "Poland"  , "North Korea"  , "Madeira Islands"  , "Azores"  , "Portugal"  , "Palestine"  , "French Polynesia"  , "Qatar"  , "Romania"  , "Russia"  , "Rwanda"  , "Western Sahara"  , "Saudi Arabia"  , "Sudan"  , "South Sudan"  , "Senegal"  , "Singapore"  , "South Sandwich Islands"  , "South Georgia"  , "Saint Helena"  , "Ascension Island"  , "Solomon Islands"  , "Sierra Leone"  , "El Salvador"  , "Somalia"  , "Serbia"  , "Slovakia"  , "Slovenia"  , "Sweden"  , "Swaziland"  , "Seychelles"  , "Syria"  , "Turks and Caicos Islands"  , "Chad"  , "Togo"  , "Thailand"  , "Tajikistan"  , "Turkmenistan"  , "Timor-Leste"  , "Tonga"  , "Trinidad"  , "Tobago"  , "Tunisia"  , "Turkey"  , "Taiwan"  , "Tanzania"  , "Uganda"  , "Ukraine"  , "Uzbekistan"  , "Vatican"  , "Vietnam"  , "Vanuatu"  , "Wallis and Futuna"  , "Samoa"  , "Yemen"  , "South Africa"  , "Zambia"  , "Zimbabwe")

map <- map_data("world", regions=regions)

#ggplot map function
plot_map <- ggplot(map, aes(x = long, y = lat, group=group)) +
    geom_polygon(fill="gray95", colour = "gray60") +
    labs(x = "", y = "") +
    theme_bw()

#define butterfly range or subsets of the range
whole_range <- coord_map("mercator", xlim = c(-30,160), ylim = c(-50,50))


Africa_range <- coord_map("mercator", xlim = c(-20,60), ylim = c(-40,40))


###############################  GBIF RECORDS  ########################################


records_GBIF <- read.csv("GBIF_Dchrysippus_20201206_20210721_20211101_20220111.csv", header = T, stringsAsFactors = F)

records_GBIF[records_GBIF == ""] <- NA

records_GBIF <- subset(records_GBIF, is.na(hindwingWhite)==FALSE | is.na(forewingTip)==FALSE | is.na(groundColour)==FALSE)

records_GBIF <- subset(records_GBIF, is.na(decimalLatitude)==FALSE & is.na(decimalLongitude)==FALSE)

table(records_GBIF$sex)
table(records_GBIF$hindwingWhite)
table(records_GBIF$forewingTip)
table(records_GBIF$groundColour)

records_GBIF$aa <- ifelse(records_GBIF$hindwingWhite == "present", 1, 0)
records_GBIF$AA <- ifelse(records_GBIF$hindwingWhite == "absent", 1, 0)
records_GBIF$Aa <- ifelse(records_GBIF$hindwingWhite == "partial", 1, 0)
records_GBIF$A. <- records_GBIF$Aa + records_GBIF$AA

records_GBIF$bb <- ifelse(records_GBIF$groundColour == "light", 1, 0)
records_GBIF$BB <- ifelse(records_GBIF$groundColour == "dark", 1, 0)
records_GBIF$Bb <- ifelse(records_GBIF$groundColour == "intermediate", 1, 0)
records_GBIF$B. <- records_GBIF$Bb + records_GBIF$BB

records_GBIF$cc <- ifelse(records_GBIF$forewingTip == "present", 1, 0)
records_GBIF$CC <- ifelse(records_GBIF$forewingTip == "absent", 1, 0)
records_GBIF$Cc <- ifelse(records_GBIF$forewingTip == "partial", 1, 0)
records_GBIF$C. <- records_GBIF$Cc + records_GBIF$CC

records_GBIF$male <- ifelse(records_GBIF$sex == "male", 1, 0)
records_GBIF$female <- ifelse(records_GBIF$sex == "female", 1, 0)


#############################  SM lab records  ##########################################

records_SM <- read.csv("SM_Dchrysippus_20231113.csv", header = T, stringsAsFactors = F)

records_SM[records_SM == ""] <- NA

records_SM$sex[records_SM$sex == "?"] <- NA

records_SM <- subset(records_SM, is.na(hindwingWhite)==FALSE | is.na(forewingTip)==FALSE | is.na(groundColour)==FALSE)

records_SM <- subset(records_SM, is.na(decimalLatitude)==FALSE & is.na(decimalLongitude)==FALSE)


table(records_SM$sex)
table(records_SM$hindwingWhite)
table(records_SM$forewingTip)
table(records_SM$groundColour)

records_SM$aa <- ifelse(records_SM$hindwingWhite == "present", 1, 0)
records_SM$AA <- ifelse(records_SM$hindwingWhite == "absent", 1, 0)
records_SM$Aa <- ifelse(records_SM$hindwingWhite == "partial", 1, 0)
records_SM$A. <- records_SM$Aa + records_SM$AA

records_SM$bb <- ifelse(records_SM$groundColour == "light", 1, 0)
records_SM$BB <- ifelse(records_SM$groundColour == "dark", 1, 0)
records_SM$Bb <- ifelse(records_SM$groundColour == "intermediate", 1, 0)
records_SM$B. <- records_SM$Bb + records_SM$BB

records_SM$cc <- ifelse(records_SM$forewingTip == "present", 1, 0)
records_SM$CC <- ifelse(records_SM$forewingTip == "absent", 1, 0)
records_SM$Cc <- ifelse(records_SM$forewingTip == "partial", 1, 0)
records_SM$C. <- records_SM$Cc + records_SM$CC

records_SM$male <- ifelse(records_SM$sex == "male", 1, 0)
records_SM$female <- ifelse(records_SM$sex == "female", 1, 0)


############################  ABRI records  ###########################################

records_ABRI <- read.csv("ABRI_Dchrysippus_selection_2019.csv", header = T, stringsAsFactors = F)

records_ABRI <- subset(records_ABRI, is.na(decimalLatitude)==FALSE & is.na(decimalLongitude)==FALSE)


table(records_ABRI$sex)
table(records_ABRI$hindwingWhite)
table(records_ABRI$forewingTip)
table(records_ABRI$groundColour)

records_ABRI$aa <- ifelse(records_ABRI$hindwingWhite == "present", 1, 0)
records_ABRI$AA <- ifelse(records_ABRI$hindwingWhite == "absent", 1, 0)
records_ABRI$Aa <- ifelse(records_ABRI$hindwingWhite == "partial", 1, 0)
records_ABRI$A. <- records_ABRI$Aa + records_ABRI$AA

records_ABRI$bb <- ifelse(records_ABRI$groundColour == "light", 1, 0)
records_ABRI$BB <- ifelse(records_ABRI$groundColour == "dark", 1, 0)
records_ABRI$Bb <- ifelse(records_ABRI$groundColour == "intermediate", 1, 0)
records_ABRI$B. <- records_ABRI$Bb + records_ABRI$BB

records_ABRI$cc <- ifelse(records_ABRI$forewingTip == "present", 1, 0)
records_ABRI$CC <- ifelse(records_ABRI$forewingTip == "absent", 1, 0)
records_ABRI$Cc <- ifelse(records_ABRI$forewingTip == "partial", 1, 0)
records_ABRI$C. <- records_ABRI$Cc + records_ABRI$CC

records_ABRI$male <- ifelse(records_ABRI$sex == "male", 1, 0)
records_ABRI$female <- ifelse(records_ABRI$sex == "female", 1, 0)


############################  Kenya collaborator records  ###########################################

records_KENYA <- read.csv("Kenya_collection_20211123.csv", header = T, stringsAsFactors = F)

records_KENYA <- subset(records_KENYA, is.na(hindwingWhite)==FALSE | is.na(forewingTip)==FALSE | is.na(groundColour)==FALSE)

records_KENYA <- subset(records_KENYA, is.na(decimalLatitude)==FALSE & is.na(decimalLongitude)==FALSE)


table(records_KENYA$sex)
table(records_KENYA$hindwingWhite)
table(records_KENYA$forewingTip)
table(records_KENYA$groundColour)

records_KENYA$aa <- ifelse(records_KENYA$hindwingWhite == "present", 1, 0)
records_KENYA$AA <- ifelse(records_KENYA$hindwingWhite == "absent", 1, 0)
records_KENYA$Aa <- ifelse(records_KENYA$hindwingWhite == "partial", 1, 0)
records_KENYA$A. <- records_KENYA$Aa + records_KENYA$AA

records_KENYA$bb <- ifelse(records_KENYA$groundColour == "light", 1, 0)
records_KENYA$BB <- ifelse(records_KENYA$groundColour == "dark", 1, 0)
records_KENYA$Bb <- ifelse(records_KENYA$groundColour == "intermediate", 1, 0)
records_KENYA$B. <- records_KENYA$Bb + records_KENYA$BB

records_KENYA$cc <- ifelse(records_KENYA$forewingTip == "present", 1, 0)
records_KENYA$CC <- ifelse(records_KENYA$forewingTip == "absent", 1, 0)
records_KENYA$Cc <- ifelse(records_KENYA$forewingTip == "partial", 1, 0)
records_KENYA$C. <- records_KENYA$Cc + records_KENYA$CC

records_KENYA$male <- ifelse(records_KENYA$sex == "male", 1, 0)
records_KENYA$female <- ifelse(records_KENYA$sex == "female", 1, 0)


##############################  David Smith Museum and Field records  #########################################

records_DAS <- read.csv("Dchrysippus_DAS_records_Auguest2021_curated.csv", header = T, stringsAsFactors = F)

#first turn all NA vlalues to zero, so we can add without problems
for (trait in c("aa", "Aa", "A.", "bb", "B.", "cc", "Cc", "C.")){
    records_DAS[is.na(records_DAS[,trait]),trait] <- 0
    }

# Make entries for known hets as well as all dominant
records_DAS$AA <- records_DAS$A.
records_DAS$A. <- records_DAS$AA + records_DAS$Aa

records_DAS$CC <- records_DAS$C.
records_DAS$C. <- records_DAS$CC + records_DAS$Cc

#the DAS records do not include any Bb, so just record that as zero
records_DAS$BB <- records_DAS$B.
records_DAS$Bb <- 0

#here we are working with individuals, so we only take the sites where a single individual is recorded

single_ind <- which(apply(records_DAS[,c("bb", "B.")], 1, sum) ==1 &
                    apply(records_DAS[,c("cc", "C.")], 1, sum) ==1 &
                    (records_DAS$sex == "m" | records_DAS$sex == "f"))

records_DAS <- records_DAS[single_ind,]

###################################  combine all records  ####################################


records <- rbind(records_GBIF[,c("decimalLongitude", "decimalLatitude","AA", "Aa", "A.", "aa", "BB", "Bb", "B.", "bb", "CC", "Cc", "C.", "cc")],
                 records_SM[,c("decimalLongitude", "decimalLatitude","AA", "Aa", "A.", "aa", "BB", "Bb", "B.", "bb", "CC", "Cc", "C.", "cc")],
                 records_ABRI[,c("decimalLongitude", "decimalLatitude","AA", "Aa", "A.", "aa", "BB", "Bb", "B.", "bb", "CC", "Cc", "C.", "cc")],
                 records_KENYA[,c("decimalLongitude", "decimalLatitude","AA", "Aa", "A.", "aa", "BB", "Bb", "B.", "bb", "CC", "Cc", "C.", "cc")],
                 records_DAS[,c("decimalLongitude", "decimalLatitude","AA", "Aa", "A.", "aa", "BB", "Bb", "B.", "bb", "CC", "Cc", "C.", "cc")])


for (A_geno in c("AA", "Aa", "aa")){
    for (B_geno in c("BB", "Bb", "bb")){
        for (C_geno in c("CC", "Cc", "cc")){
            records[,paste0(A_geno, B_geno, C_geno)] <- ifelse(records[,A_geno] == 1 & records[,B_geno] == 1 & records[,C_geno] == 1, 1, 0)
            }
        }
    }


for (A_geno in c("A.", "aa")){
    for (B_geno in c("B.", "bb")){
        for (C_geno in c("C.", "cc")){
            records[,paste0(A_geno, B_geno, C_geno)] <- ifelse(records[,A_geno] == 1 & records[,B_geno] == 1 & records[,C_geno] == 1, 1, 0)
            }
        }
    }

records$morph_definable <- ifelse(is.na(records$aa) | is.na(records$bb) | is.na(records$cc), 0, 1)


#forewing only
for (B_geno in c("BB", "Bb", "bb")){
    for (C_geno in c("CC", "Cc", "cc")){
        records[,paste0(B_geno, C_geno)] <- ifelse(records[,B_geno] == 1 & records[,C_geno] == 1, 1, 0)
        }
    }

for (B_geno in c("B.", "bb")){
    for (C_geno in c("C.", "cc")){
        records[,paste0(B_geno, C_geno)] <- ifelse(records[,B_geno] == 1 & records[,C_geno] == 1, 1, 0)
        }
    }

records$FW_definable <- ifelse(is.na(records$bb) | is.na(records$cc), 0, 1)


records$morph_dom <- ifelse(records$aabbcc == 1, "alcippus",
                        ifelse(records$A.B.cc == 1, "orientis",
                               ifelse(records$A.bbC. == 1, "klugii",
                                      ifelse(records$A.bbcc == 1, "chrysippus",
                                             ifelse(records$morph_definable == 1, "intermediate", NA)))))

records$morph <- ifelse(records$aabbcc == 1, "alcippus",
                        ifelse(records$AABBcc == 1, "orientis",
                            ifelse(records$AAbbCC == 1, "klugii",
                                    ifelse(records$AAbbcc == 1, "chrysippus",
                                            ifelse(records$morph_definable == 1, "intermediate", NA)))))

records$FW_dom <- ifelse(records$bbcc == 1, "chrysippus",
                        ifelse(records$B.cc == 1, "orientis",
                               ifelse(records$bbC. == 1, "klugii",
                                      ifelse(records$FW_definable == 1, "intermediate", NA))))


records$FW <- ifelse(records$bbcc == 1, "chrysippus",
                        ifelse(records$BBcc == 1, "orientis",
                               ifelse(records$bbCC == 1, "klugii",
                                      ifelse(records$FW_definable == 1, "intermediate", NA))))


morph_cols <- c(alcippus="#5db463", orientis="#2143d1", klugii="#ffac07", chrysippus="#bc4754" ,intermediate="black")


###########################  block data  #############################

blocksize = 5

#first we define all possible blocks in the world, and get just their latitude and longitude 'name' and coordinates
block_data <- define_all_blocks(blocksize=blocksize)

#next, for each block, we identify the records that fall within that block
rows_by_block <- group_records_into_blocks(records$decimalLongitude, records$decimalLatitude,
                                           blocksize=blocksize, block_data=block_data)


#for each trait below, the records file has a count of number of individuals with that trait
# We have already recorded which records fall with in each block, so now we just sum
# the number of individuals with each trait for each block

traits = c("bbcc", "BBcc", "bbCC", "FW_definable")

counts_by_block <- t(sapply(rows_by_block, function(rows) apply(records[rows, traits], 2, sum, na.rm=T)))

#add these counts to the block data frame
block_data <- cbind(block_data, counts_by_block)

#frequencies of FW morphs

morphs <- c("bbcc", "BBcc", "bbCC")
for (morph in morphs){
    block_data[,paste0("freq_FW_",morph)] <- block_data[,morph] / block_data$FW_definable
    }

block_data$freq_FW_other <- 1 - apply(block_data[,paste0("freq_FW_",morphs)], 1, sum)


############################# plot map FW morphs ######################################

block_data_subset <- subset(block_data, FW_definable>0)

columns = paste0("freq_FW_", c("bbcc", "BBcc", "bbCC", "other"))


# png("morph_pie_map.png", width=3000, height=2000, res=300, bg=NA)
svg("FW_pie_map.svg", width=10, height=5, bg=NA)
plot_map + whole_range +
    geom_scatterpie(aes(x = block_midpoint_long, y = block_midpoint_lat, r=cap(FW_definable,5)*.5),
                    data = block_data_subset, cols = columns, color=NA, show.legend=FALSE) +
    scale_fill_manual(values=c("#bc4754","#2143d1","#ffac07" ,"black"))
dev.off()



#########################################################################################
#################################   points   ############################################
#########################################################################################


############# plot forewing points jittered ##################

svg("FW_map.svg", width=10, height=5, bg=NA)
plot_map + whole_range +
    geom_jitter(aes(x = decimalLongitude, y = decimalLatitude, colour=FW),
                data = records[is.na(records$FW) == FALSE,],
                inherit.aes=FALSE, width=0.2, height=0.2, alpha=1, size=1, shape=16, show.legend=FALSE) + 
    scale_color_manual(values = morph_cols)
dev.off()



############################### MAP DC174 sequenced samples ##############################################

morph_cols <- c(alcippus="#5db463", orientis="#2143d1", klugii="#ffac07", chrysippus="#bc4754", intermediate="black", unknown="gray")


DC174 <- read.csv("DC174_data.csv", as.is=T)

DC174$FW <- ifelse(DC174$groundColour == "light" & DC174$forewingTip == "present", "chrysippus",
                   ifelse(DC174$groundColour == "light" & DC174$forewingTip == "absent", "klugii",
                          ifelse(DC174$groundColour == "dark" & DC174$forewingTip == "present", "orientis",
                                 ifelse(DC174$groundColour == "intermediate" | DC174$forewingTip == "partial", "intermediate", "unknown"))))

DC174_sub <- subset(DC174, is.na(lat) == FALSE & is.na(long) == FALSE & is.na(FW) == FALSE)

svg("DC174_map.svg", height=5, bg=NA)
plot_map + Africa_range +
    geom_jitter(aes(x = long, y = lat, colour=FW, shape=group),
                data = DC174_sub,
                inherit.aes=FALSE, width=.5, height=.5, alpha=1, size=1.5, stroke=.7, show.legend=FALSE) + 
    scale_color_manual(values = morph_cols) + 
    scale_shape_manual(values = c(4, 4))

dev.off()


#### DC174 locations without jitter

DC174 <- read.csv("DC174_data.csv", as.is=T)

DC174_sub <- subset(DC174, is.na(lat) == FALSE & is.na(long) == FALSE)

unique_coords <- unique(apply(DC174_sub[,c("lat","long")], 1, paste, collapse="_"))

unique_coords_df <- as.data.frame(t(simplify2array(strsplit(unique_coords, "_"))), stringsAsFactors=FALSE)

names(unique_coords_df) <- c("lat", "long")

plot_map + Africa_range +
    geom_point(aes(x = long, y = lat),
                data = DC174_sub,
                inherit.aes=FALSE, size=1, stroke=.7, show.legend=FALSE)

ggsave("distributions/DC174_Africa.pdf", width=3, height=4, bg='transparent')
