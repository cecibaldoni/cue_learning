#script for spatial analysis of shrew movement
#Elham Nourani, PhD. May 19.2022. Konstanz, Germany.
#update to analyse all trials in all seasons: July 4. 2022
#update Oct. 12. 2022
#--------------------------------------------------------------------------

library(tidyverse)
library(sf) #for spatial manipulation
library(mapview) #for interactive visualization of the spatial data
library(parallel)
library(ggplot2)
# STEP 1: open data for all trials ------------------------------------------

#doors <- read.csv("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/jul4_data/trial_door.csv") %>% 
#  mutate(Trial = paste0("T",trial_n))
doors <- read.csv("~/Documents/data/cue_learning/trial_door.csv") %>% 
  mutate(Trial = paste0("T",trial_n))
# coords <- read.csv("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/oct12_data/Food_door coordinates.csv") %>% 
#   mutate(Trial = paste0("T",TRIAL)) %>% 
#   mutate(unique_trial_ID = paste(SEASON, Trial, ID, sep = "_"))

coords <- read.csv("~/Documents/data/cue_learning/food_door_coordinates.csv", sep = ";",header = TRUE) %>%
  mutate(Trial = paste0("T",TRIAL)) %>%
  mutate(unique_trial_ID = paste(SEASON, Trial, ID, sep = "_"))

# tracking <- lapply(list.files("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/oct12_data/tracking", full.names = T), read.csv) %>% 
#   reduce(rbind) %>% 
#   drop_na(frame) %>% 
#   mutate(unique_trial_ID = paste(Season, Trial, ID, sep = "_")) #create unique trial IDs. there are 10 trials per season per shrew

tracking <- lapply(list.files("~/Documents/data/cue_learning/tracking", full.names = T), read.csv) %>%
  reduce(rbind) %>%
  drop_na(frame) %>%
  mutate(unique_trial_ID = paste(Season, Trial, ID, sep = "_")) #create unique trial IDs. there are 10 trials per season per shrew

trialIDs=unique(tracking$unique_trial_ID)
#in tracking.csv but not in food_door_coordinates.csv
trialIDs[!(trialIDs %in% coords$unique_trial_ID)]
#in food_door_coordinates.csv but not in tracking.csv  
coords$unique_trial_ID[!(coords$unique_trial_ID %in% trialIDs)]

# # STEP 2: create spatial objects ------------------------------------------
# 
# #spatial point for food location
# food_pts <- st_as_sf(coords, coords =  c("FOOD_x","FOOD_y"))
# # # 
# # # #create a buffer of 6 cm around the food
# food_buff <- food_pts %>% 
#   st_buffer(dist = 6) 
# 
# # # #doors list: points and buffered points
# doors_pts <- coords %>%
#   st_as_sf(coords = c("A_x", "A_y")) %>%
#   st_as_sf(coords = c("B_x", "B_y")) %>%
#   st_as_sf(coords = c("C_x", "C_y")) %>%
#   st_as_sf(coords = c("D_x", "D_y")) %>%
#   mutate(shortest_path = st_distance(x = .,y = food_pts)) #length of the shortest trajectory from each door to food))
# doors_pts <- coords %>%
#   st_as_sf(coords = c("A_x", "A_y")) %>%
#   st_as_sf(coords = c("B_x", "B_y")) %>%
#   st_as_sf(coords = c("C_x", "C_y")) %>%
#   st_as_sf(coords = c("D_x", "D_y")) %>%
#   mutate(shortest_path = st_distance(x = .,y = food_coords))
# 
# doors_buff <- doors_pts %>% 
#   st_buffer(dist = 6) #adjust this buffer as you see fit!
# 
# trials_sf <- tracking %>% 
#   drop_na(x) %>% 
#   st_as_sf(coords = c("x", "y"))

# STEP 3: spatial investigations! ------------------------------------------

#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL) #we can also add the length of time spent at the door if needed

trial_ls <- split(tracking, tracking$unique_trial_ID)

#prepare cluster for parallel computation
#mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)
mycl <- makeCluster(3)

#clusterExport(mycl, c("trial_ls", "other_door_visits", "doors_pts", "food_buff","doors_buff", "doors")) #define the variable that will be used within the ParLapply call

clusterExport(mycl, c("coords", "trial_ls", "doors")) #define the variable that will be used within the ParLapply call
 
clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
   library(sf)
   library(sp)
   library(tidyverse)
 })


b <- Sys.time()
#sp_prep <- lapply(trial_ls, function(x){ 

sp_prep <-parLapply(mycl, trial_ls, function(x){ 
  
  #to call all blocks one at a time
  #x = trial_ls[[1]]
  
  #extract food coordinates for this trial AND convert to a sf object
  food_coords <- coords %>% 
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>% 
    dplyr::select(c("FOOD_x", "FOOD_y")) %>% 
    st_as_sf(coords = c("FOOD_x", "FOOD_y"))
    
  food_buffer <- food_coords %>% 
    st_buffer(dist = 4) #half of the length of the largest possible shrew
    
  #extract door coordinates for this trial AND convert to a sf object
  trial_door_ID <- doors %>% 
    filter(Trial == unique(x$Trial)) %>% 
    pull(door)
  
  trial_door_coords <- coords %>% 
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>% 
    dplyr::select(starts_with(trial_door_ID)) %>% 
    st_as_sf(coords = colnames(.))
  
  trial_door_buffer <- trial_door_coords %>% 
    st_buffer(dist = 3.5)
  
  #convert the track into an sf object
  track_sf <- x %>% 
    st_as_sf(coords = c("x", "y"))
  
  #find the points of overlap between track and food
  #WARNING MESSAGE
  at_food <- track_sf %>% 
    st_intersection(food_buffer) %>% 
    as.data.frame() %>% 
    arrange(frame) %>% #arrange by time/frame
    mutate(timediff = frame - lag(frame)) %>% 
    mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1,0)) %>% 
    mutate(visit_seq = cumsum(new_timediff))
  
    #plot to check
  # mapview(track_sf)  + 
  #   mapview(food_buffer, color = "yellow") + 
  #   mapview(at_food, color = "yellow") + 
  #   mapview(food_coords, color = "black")
  
  #add at food journey info to x
  track_sf_2 <- track_sf %>% 
    full_join(at_food[c("frame", "visit_seq")]) %>% 
    arrange(frame) %>% 
    mutate(food_journey = ifelse(frame == head(at_food$frame,1), "arrival", #the first point in the at_food data is the arrival
                                 ifelse(frame == tail(at_food[at_food$visit_seq == 1, "frame"],1), "departure", #the last point in the first visit to food
                                        ifelse(frame < head(at_food$frame,1), "trip_to", #points before point of arrival are trip to food
                                               ifelse(between(frame, head(at_food[at_food$visit_seq == 1, "frame"],1), tail(at_food[at_food$visit_seq == 1,"frame"],1)), "at_food", #time spent at food, between arrival and departure from food
                                                      ifelse(frame %in% at_food[at_food$visit_seq != 1,"frame"], 
                                                             paste("trip_back_revisit", visit_seq, sep = "_"),
                                                             "trip_back")))))) #all other points are trip back from food
  
  #plot to check
  #mapview(track_sf_2, zcol = "food_journey") + mapview(food_buffer, color = "orange") 
  
  # check for overlap between the return trip and other doors
  #PROBLEM WITH THE CODE. "door" not found, not sure how to fix
  # other_doors <- track_sf_2 %>%
  #    filter(food_journey == "trip_back") %>%
  #    st_intersection(trial_door_buffer %>% filter(door != unique(trial_door$door))) #exclude the trial door
  # # 
  # #if there was a visit to another door, save info to other_door_visits dataframe
  # if(nrow(other_doors) > 0){
  #   new_visits <- other_doors %>%
  #     group_by(door) %>%  slice(1) %>%  #this will give you one row per other door visited
  #     dplyr::select(c("ID", "Trial", "door")) %>%
  #     st_drop_geometry() #convert back to non-spatial object
  # 
  #   #append to other_door_visits
  #   other_door_visits <<- rbind(other_door_visits,new_visits) #double arrow assignment operator allows to modify the dataframe in the global environment
  # }
  
  return(x)
})

Sys.time() - b 
stopCluster(mycl) 

mapview(trial_door_coords) + mapview(food_buffer) + mapview(food_coords) + mapview(trial_door_buffer) + mapview(track_sf)


#saveRDS(sp_prep, file = "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_sp.rds")
saveRDS(sp_prep, file = "/home/ceci/Documents/data/trials_sp.rds")
#saveRDS(other_door_visits, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_other_doors.rds")
saveRDS(other_door_visits, "/home/ceci/Documents/data/trials_other_doors.rds")

# STEP 4: Calculate speed, distance, tortuosity ------------------------------------------

#sp_prep <- readRDS("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_sp.rds")

trial_summaries <- sp_prep %>%
  reduce(rbind) %>% #merge the list elements
  group_by(unique_trial_id, food_journey) %>%  #calculate length and duration for each leg of the food journey within each trial for each individual
  arrange(frame) %>% 
  summarise(do_union = FALSE,
            duration = tail(time,1) - head(time,1)) %>%  #calculate the number of seconds spent in each part of the food journey (the difference between the last timestamp and the first)
  st_cast("LINESTRING") %>%  #convert to line object to calculate length.
  ungroup() %>% 
  mutate(length = st_length(.)) %>%  #in cm. make sure it makes sense.
  rowwise() %>% 
  mutate(tortuosity = ifelse(food_journey %in% c("trip_to","trip_back"), 
                             length/as.numeric(doors_pt[doors_pt$Door == Door,"shortest_path"])[1], NA)) %>% #calculate tortuosity as the ratio between length of track and the shortest path (only for way to food and back)
  st_drop_geometry() #remove the spatial portion of the dataset

saveRDS(trial_summaries, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_summaries.rds")

# STEP 5: Check for visit to previous trial's door ------------------------------------------

other_door_visits <- readRDS("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_other_doors.rds")

#I'm adding this info to the dataframe that associates the trial to its door
trial_door <- trial_door %>% 
  rowwise() %>% 
  mutate(visit_to_prev_door = ifelse(Trial == "T1", NA, #assign NA for the first trial
                                     ifelse(nrow(other_door_visits[other_door_visits$Trial == Trial,]) == 0, "No", #if this trial isn't in the other_door_visitis, assign no
                                     ifelse(lag(Door,1) %in% as.character(other_door_visits[other_door_visits$Trial == Trial, "Door.1"]), "Yes", "No"))))




