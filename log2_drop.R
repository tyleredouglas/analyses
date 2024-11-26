#takes hits file (local maxima) and sites file (all sites from model output) from get_loci.R
# and outputs min and max coordinates of peak area define by a 2 log drop

log2_drop <- function(hits_file, sites_file, output_name, dest) {
    
    hits <- read.table(hits_file, header = TRUE)   
    sites <- read.table(sites_file, header = TRUE)

    #initialize dataframe for results
    peak_coords <- tibble(hit_position = character(), hit_logp = numeric(), max_pos = character(), max_p = numeric(), min_pos = character(), min_p = numeric()) 
    

    for (i in 1:nrow(hits)) {
      hit_row <- hits[i,]
      hit_pos <- as.character(hit_row[1])
      hit_p <- as.numeric(hit_row[6])
      cpos <- as.character(hit_row[1])
      logp_max <- as.numeric(hit_row[6])
      logp_min <- as.numeric(hit_row[6])

      #define minimum logp to establish window boundaries (log2drop)
      drop2 <- as.numeric(logp_max - 2)
      
      #get row ID of hit in sites table
      hit_id_max <- as.numeric(sites$row_id[match(cpos, sites$chr_pos)])
      hit_id_min <- as.numeric(sites$row_id[match(cpos, sites$chr_pos)])
      
      #row row ID of the two sites before and after hit
      #need to assess sites three at a time b/c certain sites are outliers with low logp, will prematurely truncate window
      id_plus_one <- hit_id_max + 1
      id_plus_two <- hit_id_max + 2
      id_minus_one <- hit_id_min - 1
      id_minus_two <- hit_id_min - 2

      #create variables for logp of two flanking positions, temporarily defining with value of logp of hit
      logp_plus_one <- logp_max
      logp_plus_two <- logp_max
      logp_minus_one <- logp_min
      logp_minus_two <- logp_min

      #minimum logp will depend on power and effect size in given dataset, for nw_n using logp = 1.3 (p = 0.05)
      if (drop2 > 1.3) {
        
        while ((logp_max > drop2 | logp_plus_one > drop2 | logp_plus_two > drop2) & id_plus_two < nrow(sites)) {
          
          #get logp of flanking sites
          current_row <- subset(sites, row_id == hit_id_max)
      	  row_plus_one <- subset(sites, row_id == id_plus_one)
	        row_plus_two <- subset(sites, row_id == id_plus_two)
          coord_max <- as.character(current_row[1])
          logp_max <- as.numeric(current_row[4])
          logp_plus_one <- as.numeric(row_plus_one[4])
	        logp_plus_two <- as.numeric(row_plus_two[4])

          if (logp_max > 1.3 | logp_plus_one > drop2 | logp_plus_two > drop2) {
            
            #coordinate and logp of site to be recorded into output table
            max_p_rec <- logp_max
            max_coord_rec <- coord_max
          }

          #advance one position
          hit_id_max <- hit_id_max + 1
          id_plus_one <- id_plus_one + 1
      	  id_plus_two <- id_plus_two + 1

        }
        
        #same as above but in reverse direction
        while ((logp_min > drop2 | logp_minus_one > drop2 | logp_minus_two > drop2) & id_minus_two > 0)  {
          current_row <- subset(sites, row_id == hit_id_min)
          row_minus_one <- subset(sites, row_id == id_minus_one)
	        row_minus_two <- subset(sites, row_id == id_minus_two)
          coord_min <- as.character(current_row[1])
          logp_min <- as.numeric(current_row[4])
          logp_minus_one <- as.numeric(row_minus_one[4])
	        logp_minus_two <- as.numeric(row_minus_two[4])

          if (logp_min > 1.3 | logp_minus_one > drop2 | logp_minus_two > drop2) {
            min_p_rec <- logp_min
            min_coord_rec <- coord_min
          }
          hit_id_min <- hit_id_min - 1 
	        id_minus_one <- id_minus_one - 1
	        id_minus_two <- id_minus_two -1
        }
      }
      
      else {
        
        while (logp_max > 1.3) {
          current_row <- subset(sites, row_id == hit_id_max)
          coord_max <- as.character(current_row[1])
          logp_max <- as.numeric(current_row[4])
          
          if (logp_max > 1.3) {
            max_p_rec <- logp_max
            max_coord_rec <- coord_max
          }
          hit_id_max <- hit_id_max + 1
          
        }
        
        while (logp_min > 1.3) {
          current_row <- subset(sites, row_id == hit_id_min)
          coord_min <- as.character(current_row[1])
          logp_min <- as.numeric(current_row[4])
          
          if (logp_min > 1.3) {
            min_p_rec <- logp_min
            min_coord_rec <- coord_min
          }
          hit_id_min <- hit_id_min - 1 
        }
      }
      
      #store hit, hit coordinate, as well as min and max coordinate and logp of peak defined by a log2 drop
      peak_coords <- add_row(peak_coords, hit_position = hit_pos, hit_logp = hit_p,
                           max_pos = max_coord_rec, max_p = max_p_rec, 
                           min_pos = min_coord_rec, min_p = min_p_rec)
    }
    
    output_file <- peak_coords %>%
      separate(max_pos, c("max_chr", "max_coord")) %>%
      separate(min_pos, c("min_chr", "min_coord")) %>%
      select(hit_position, hit_logp, max_chr, min_coord, max_coord, min_p, max_p) 
    
    output_file$max_chr <- recode(output_file$max_chr, 
                                  "chr2L" = "2L",
                                  "chr2R" = "2R",
                                  "chr3L" = "3L",
                                  "chr3R" = "3R",
                                  "chrX" = "X")
    
    filename=paste(dest, "/", output_name, sep='')
    write.table(output_file, filename,row.names=FALSE,sep="\t",quote=FALSE)

}
