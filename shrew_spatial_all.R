#script for spatial analysis of shrew movement
#Elham Nourani, PhD. May 19.2022. Konstanz, Germany.
#update to analyse all trials in all seasons: July 4. 2022
#update Oct. 12. 2022
#--------------------------------------------------------------------------

library(tidyverse)
library(sf) #for spatial manipulation
library(mapview) #for interactive visualization of the spatial data
library(parallel)

# STEP 1: open data for all trials ------------------------------------------
#coords <- read.csv("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/jul4_data/coordinates_cue.csv") %>% 
#  mutate(season = str_replace_all(season, " ", ""))
doors <- read.csv("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/jul4_data/trial_door.csv") %>% 
  mutate(Trial = paste0("T",trial_n))

coords <- read.csv("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/oct12/Food_door coordinates.csv")

tracking <- lapply(list.files("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/oct12_data/tracking", full.names = T), read.csv) %>% 
  reduce(rbind) %>% 
  drop_na(frame)

# STEP 1: open data for 5 trials ------------------------------------------
# files_ls <- list.files("/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/may16_data/SHREWCECI/", pattern = ".csv", full.names = T)
# 
# trials <- lapply(files_ls, read.csv) %>% 
#   reduce(full_join, by = names(lapply(files_ls,read.csv)[[1]])) %>% 
#   dplyr::select(c(1:11)) %>%  #the 5th trial has two extra columns. removing these columns from the final dataset
#   drop_na(x) #remove rows with NAs. maybe the merging procedure wasn't super efficient
# 
# food_location <- data.frame(x = 127.63317, y = 26.34548)
# 
# doors <- data.frame(Door = c("DoorA", "DoorB", "DoorC", "DoorD"),
#                     x = c(50.19882, 103.95784, 159.67497, 106.44998),
#                     y = c(75.12022, 18.33503, 72.62808, 127.98919))
# 
# trial_door <- data.frame(Trial = c("T1", "T2", "T3", "T4", "T5"),
#                          Door = c("DoorA", "DoorC", "DoorA", "DoorB", "DoorD"))
# 
# #add trial door to dataset
# trials <- trials %>% 
#   full_join(trial_door, by = "Trial")

# STEP 2: create spatial objects ------------------------------------------

#spatial point for food location
food_pts <- st_as_sf(coords, coords =  c("food_x","food_y"))

#create a buffer of 6 cm around the food
food_buff <- food_pts %>% 
  st_buffer(dist = 6) 

#doors list: points and buffered points
doors_pts <- coords %>% 
  st_as_sf(coords = c("x", "y"))# %>% 
  #mutate(shortest_path = st_distance(x = .,y = food_pt)) #length of the shortest trajectory from each door to food))

doors_buff <- doors_pts %>% 
  st_buffer(dist = 6) #adjust this buffer as you see fit!

trials_sf <- tracking %>% 
  drop_na(x) %>% 
  st_as_sf(coords = c("x", "y"))

# STEP 3: spatial investigations! ------------------------------------------

#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL) #we can also add the length of time spent at the door if needed

trials_sf$unique_trial_id <- paste(trials_sf$Season, trials_sf$Trial, trials_sf$ID, sep = "_")

trial_ls <- split(trials_sf, trials_sf$unique_trial_id)

#prepare cluster for parallel computation
mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)

clusterExport(mycl, c("trial_ls", "other_door_visits", "doors_pts", "food_buff","doors_buff", "doors")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
})


b <- Sys.time()
#sp_prep <- lapply(trial_ls, function(x){ 
sp_prep <-parLapply(mycl, trial_ls, function(x){ 
  
  #extract door info for this trial
  trial_door <- doors_pts %>% 
    filter(ID == unique(x$ID) & season == unique(x$Season) & door == doors[doors$Trial == unique(x$Trial), "door"]) 
  
  trial_food <- food_buff[food_buff$ID == unique(x$ID) & food_buff$season == unique(x$Season) & food_buff$door == trial_door$door,]
  
  #extract location of all the doors for trial x
  trial_door_buff <- doors_buff[doors_buff$ID == unique(x$ID) & doors_buff$season == unique(x$Season),]
  
  #plot track and door and food and other doors!
  #mapview(x) + mapview(trial_door, color = "red") + mapview(trial_food, color = "orange") + mapview(trial_door_buff, color = "black")
  
  #find the points of overlap between track and food
  at_food <- x %>% 
    st_intersection(trial_food) %>% 
    arrange(frame) #arrange by time/frame
  
  #plot to check
  #mapview(x)  + mapview(trial_food, color = "orange") + mapview(at_food, color = "yellow")
  
  #add at food journey info to x
  x <- x %>% 
    mutate(food_journey = ifelse(frame == head(at_food$frame,1), "arrival", #the first point in the at_food data is the arrival
                                 ifelse(frame == tail(at_food$frame,1), "departure", #the last point in the at_food data is the arrival
                                        ifelse(frame < head(at_food$frame,1), "trip_to", #points before point of arrival are trip to food
                                               ifelse(between(frame, head(at_food$frame,1), tail(at_food$frame,1)), "at_food", #time spent at food, between arrival and departure from food
                                                      "trip_back"))))) #all other points are trip back from food
  
  #plot to check
  #mapview(x, zcol = "food_journey") + mapview(trial_food, color = "orange") 
  
  #check for overlap between the return trip and other doors
  other_doors <- x %>% 
    filter(food_journey == "trip_back") %>% 
    st_intersection(trial_door_buff %>%  filter(door != unique(trial_door$door))) #exclude the trial door
  
  #if there was a visit to another door, save info to other_door_visits dataframe
  if(nrow(other_doors) > 0){
    new_visits <- other_doors %>% 
      group_by(door) %>%  slice(1) %>%  #this will give you one row per other door visited
      dplyr::select(c("ID", "Trial", "door")) %>% 
      st_drop_geometry() #convert back to non-spatial object
    
    #append to other_door_visits
    other_door_visits <<- rbind(other_door_visits,new_visits) #double arrow assignment operator allows to modify the dataframe in the global environment
  }
  
  return(x)
})

Sys.time() - b 
stopCluster(mycl) 


#saveRDS(sp_prep, file = "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_sp.rds")
#saveRDS(other_door_visits, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_other_doors.rds")


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
  mutate(visit_to_prev_door = ifelse(Trial == "T1", NA, #assign AN for the first trial
                                     ifelse(nrow(other_door_visits[other_door_visits$Trial == Trial,]) == 0, "No", #if this trial isn't in the other_door_visitis, assign no
                                     ifelse(lag(Door,1) %in% as.character(other_door_visits[other_door_visits$Trial == Trial, "Door.1"]), "Yes", "No"))))




