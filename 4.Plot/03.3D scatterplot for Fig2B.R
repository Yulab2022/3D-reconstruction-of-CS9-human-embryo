
library(rgl)
#install.packages("readobj")
library(readobj) 
library(RColorBrewer)
library(plotly)
#install.packages("geomorph")
library(geomorph)
library(Seurat)


data.shell=function(obj){
  shell =  list(vt = data.frame(x = obj[[1]]$positions[3,],
                                y = obj[[1]]$positions[2,],
                                z = obj[[1]]$positions[1,]),
                i = obj[[1]]$indices[1,], 
                j = obj[[1]]$indices[2,], 
                k = obj[[1]]$indices[3,])
  return(shell)
}

##import data
all_nt =data.shell(read.obj("NT_module.obj")[[1]]) ##3D module of neural tube
all_somite =data.shell(read.obj("Somite_module.obj")[[1]]) ##3D module of somite
data3_nt_meta <- read.csv("NT_paste.csv") ##3D point cloud of neural tube


p1 <- plot_ly(z = data3_nt_meta$X_xz_02_fz_90/100, y = -(data3_nt_meta$Y_xz_02_90/100), x = data3_nt_meta$Z_paste * 0.8, 
              type = "scatter3d", mode = "markers", color=factor(data3_nt_meta$NT_celltype),
              colors = class_color_map3,
              marker = list(size=2.5, line = list(opacity = 0.8)),
              alpha = 1)%>% 
  add_trace(type = 'mesh3d',
            data = all_nt$vt,
            x = all_nt$vt$y,
            y = all_nt$vt$x,
            z = all_nt$vt$z,
            i = all_nt$i, j = all_nt$j, k = all_nt$k,
            facecolor = rep("gray",length(all_nt$i)),
            opacity = 0.015,
            name = "shell"
  ) %>%
  add_trace(type = 'mesh3d',
            data = all_somite$vt,
            x = all_somite$vt$y,
            y = all_somite$vt$x,
            z = all_somite$vt$z,
            i = all_somite$i, j = all_somite$j, k = all_somite$k,
            facecolor = rep("#BDA7CB",length(all_somite$i)),
            opacity = 0.015,
            name = "shell"
  ) %>%
  layout(scene=list(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,title=""),
                    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,title=""),
                    zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,title="")),
         legend = list(font = list(size = 15),itemwidth = 20,itemsizing = "constant",
                       itemlist = list(marker = list(size = 40))))

p1
