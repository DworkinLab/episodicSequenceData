source('Pi_PlotFunction.R')

# Example for one file (can repeat and change as necessary)
# Input is file: and the mapper for the title of the plot
Pi_PlotFunction('F115ConR1_TAGCTT_novo.pi', "Novoalign")


# Example to run all in a for loop
# Creates a list of all .pi files (from mapper)
MyPi <- list.files(pattern="novo.pi")
Mapper <- "Novoalign"

for (file in MyPi){
  Pi_PlotFunction(file, Mapper)
}
