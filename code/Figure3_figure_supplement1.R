```{r}
#read the saved visNetwork file from source_data/ ------------------------------

conn_grouped_graph.visn <- readRDS("source_data/Figure3_figure_supplement1_source_data1.rds")
```

```{r}
# plot grouped graph ----------------------------------

#overwrite group value (partition) with side value or other value (for colouring)
conn_grouped_graph.visn$nodes$group <-
  as.character(conn_grouped_graph.visn$nodes$class)

## copy column "weight" to new column "value" in list "edges"
conn_grouped_graph.visn$edges$value <- 
  sqrt(conn_grouped_graph.visn$edges$weight)

## copy column "weight" to new column "value" in list "nodes"
conn_grouped_graph.visn$nodes$value <- 
  conn_grouped_graph.visn$nodes$weight

#for plotting with different colors, remove colour info (which takes precedence over group colour)
conn_grouped_graph.visn$nodes$color <- c()
Reds <- brewer.pal(9, 'Reds')

```

```{r}
# shorten names and assign names to node labels --------------------------------

#shorten names after _ and \s
conn_grouped_graph.visn$nodes$CATMAID_name  <- sub("_.*$", "", conn_grouped_graph.visn$nodes$CATMAID_name)
conn_grouped_graph.visn$nodes$CATMAID_name  <- sub("\\s.*$", "", conn_grouped_graph.visn$nodes$CATMAID_name)

#assign name to label
conn_grouped_graph.visn$nodes$label <-  conn_grouped_graph.visn$nodes$CATMAID_name

```

```{r}
# generate visNet graph ---------------------------

visNet <- visNetwork(conn_grouped_graph.visn$nodes, 
                     conn_grouped_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE) %>%
  visPhysics(solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -150000,
                                       centralGravity = 10)) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=1, max=25),
             color = list(color = "#454545", opacity = 0.2)) %>%
  visNodes(borderWidth=0.3, 
             color = list(border='black'),
             scaling=list(min=20, max=50),
             font = list(size = 35, 
                         color = "black", 
                         background = "#DEDEDE")) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE)  %>%
  visGroups(groupname = "muscle", color='grey', shape = "dot", 
            opacity=0.5) %>%
  visGroups(groupname = "epithelia_cell", color="#0072B2", shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "circumchaetal", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "circumacicular", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "acicula", color="#D55E00", shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "chaeta", color="#D55E00", shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "hemichaetal", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "ciliated cell", color="#E69F00", shape = "dot", 
                                     opacity=0.8) 

visNet
```

```{r}
# save graph ---------------------------

#save as html
saveNetwork(visNet, "pictures/Grouped_desmo_connectome.html", selfcontained = TRUE)
#create png snapshot
webshot2::webshot(url="pictures/Grouped_desmo_connectome.html",
                    file="pictures/Grouped_desmo_connectome.png",
                    vwidth = 1500, vheight = 1500, #define the size of the browser window
                    cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
```



