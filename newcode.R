trialIDs=unique(temp_track$unique_trial_ID)

for(i in 1:length(trialIDs)){
  trial_needed=trialIDs[i]
  temp_track=tracking[which(tracking$unique_trial_ID==trial_needed),]
  temp_track$z=temp_track$x+1i*temp_track$y
  
  doors_temp=coords[which(coords$unique_trial_ID==trial_needed),]
  doors_temp$A_x=doors_temp$A_x*doors_temp$cmxpixel
  doors_temp$A_y=doors_temp$A_y*doors_temp$cmxpixel
  doors_temp$B_x=doors_temp$B_x*doors_temp$cmxpixel
  doors_temp$B_y=doors_temp$B_y*doors_temp$cmxpixel
  doors_temp$C_x=doors_temp$C_x*doors_temp$cmxpixel
  doors_temp$C_y=doors_temp$C_y*doors_temp$cmxpixel
  doors_temp$D_x=doors_temp$D_x*doors_temp$cmxpixel
  doors_temp$D_y=doors_temp$D_y*doors_temp$cmxpixel
  doors_temp$FOOD_x=doors_temp$FOOD_x*doors_temp$cmxpixel
  doors_temp$FOOD_y=doors_temp$FOOD_y*doors_temp$cmxpixel
  
  Food=doors_temp$FOOD_x+1i*doors_temp$FOOD_y
  temp_track$distance_to_food=temp_track$z-Food
  temp_track$distance_to_food=Mod(temp_track$distance_to_food)
  arrival=temp_track[which(temp_track$distance_to_food<=3),]
  arrival=arrival[which(arrival$frame==min(arrival$frame,na.rm=TRUE)),]
  temp_track$arrival=NA
  temp_track$arrival[which(temp_track$frame<unique(arrival$frame))]="before"
  temp_track$arrival[which(temp_track$frame>=unique(arrival$frame))]="after"
  
  ggplot()+geom_path(data=temp_track, aes(x=x, y=y, color=frame))+
    annotate(
      "text", label = "A",
      x=doors_temp$A_x, y=doors_temp$A_y,  size = 8, colour = "red"
    ) +
    annotate(
      "text", label = "B",
      x=doors_temp$B_x, y=doors_temp$B_y, size = 8, colour = "red"
    ) +
    annotate(
      "text", label = "C",
      x=doors_temp$C_x, y=doors_temp$C_y, size = 8, colour = "red"
    ) +
    annotate(
      "text", label = "D",
      x= doors_temp$D_x, y=doors_temp$D_y, size = 8, colour = "red"
    ) +
    geom_point(data=doors_temp, aes(x=FOOD_x, y=FOOD_y, size=2), color=I("orange"))+
    theme_classic()+ coord_equal()
  
  meanstep=mean(Mod(diff(z)),na.rm=TRUE)
}
