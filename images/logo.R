install.packages("hexSticker")

library(hexSticker)
library(magick)

imgurl <- "C:/Users/bald_local/Desktop/sdmPermormance_temp/logo_rpackage.png"
# Read the image with magick
logo <- image_read(imgurl)
s=sticker(subplot=logo,
        package="poEvaluationMetrics",
          p_size = 12,
        p_y=1.49,
          s_x = 1,       # horizontal shift
          s_y = 0.9,    # vertical shift
          s_width = 1, # size of subplot
        s_height= 1,
        h_fill= "white",
        h_color= "#38b6ff",
        p_color = "#800000",
        filename ="C:/Users/bald_local/Desktop/sdmPermormance_temp/logo.png"        )

s
