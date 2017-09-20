## rainbow colours by ggthemr
# load ggplot as ggthemr uses it's colours
library(ggplot2)
library(ggthemr)

## If dont have ggthemr already, can install ggthemr using this;
## devtools::install_github('cttobin/ggthemr')

# 18 colours about 

dark_cols <- c("gray22", "chocolate4", "slateblue", "olivedrab")

DarkCols1 <- c("#555555", dark_cols)
# remove previous effects:
ggthemr_reset()
# Define colours for your figures with define_palette
darkCols <- define_palette(
  swatch = DarkCols1, # colours for plotting points and bars
  gradient = c(lower = DarkCols1[1L], upper = DarkCols1[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)
# set the theme for your figures:
ggthemr(darkCols)


