# =======================================================================================
# Bash script to convert a set of frames stored in .png to video animation in .mp4
# Author: Lucas Berg
# =======================================================================================
#!/bin/bash

# Variables
FILENAME="frames/frame"
FRAME_RATE="250"
END_FRAME="2500"
OUTPUT_VIDEO_FILENAME="videos/elizabeth_dla_guided"
RESOLUTION="1980x1024"

# Execute the converting command using FFMPEG
ffmpeg -r ${FRAME_RATE} -f image2 -s ${RESOLUTION} -start_number 1 -i ${FILENAME}.%04d.png -vframes ${END_FRAME} -vcodec libx264 -crf 25  -pix_fmt yuv420p ${OUTPUT_VIDEO_FILENAME}.mp4

# Working version for sending .mp4 via WhatsApp
#ffmpeg -i ${OUTPUT_VIDEO_FILENAME}.mp4 -c:v libx264 -b:v 1500k -c:a aac fixedvideo.mp4
