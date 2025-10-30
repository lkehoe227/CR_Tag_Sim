library(leaflet)
library(sf)
library(dplyr)
library(tidyr)
library(tigris)
library(leaflet.extras2)
library(mapview)

options(tigris_use_cache = TRUE)

flc <- counties(state = "FL", cb = TRUE, class = "sf") |>
  st_transform(4326)

region_defs <- list(
  "MRIP 1" = c("Escambia","Santa Rosa","Okaloosa","Walton","Bay","Gulf"),
  "MRIP 2" = c("Hernando","Pasco","Wakulla","Jefferson","Taylor","Dixie","Levy","Citrus","Franklin"),
  "MRIP 3" = c("Sarasota","Pinellas","Hillsborough","Manatee"),
  "MRIP 4" = c("Charlotte","Lee","Collier"),
  "MRIP 5" = c("Miami-Dade","Broward","Palm Beach"),
  "MRIP 6" = c("St. Lucie","Indian River","Brevard","Martin"),
  "MRIP 7" = c("Flagler","St. Johns","Duval","Nassau","Volusia"),
  "MRIP 8" = c("Monroe")   # Keys
)

region_lookup <- tibble(region = names(region_defs), counties = region_defs) |>
  tidyr::unnest_longer(counties, values_to = "NAME") |>
  dplyr::select(NAME, region)

flc <- flc |> left_join(region_lookup, by = "NAME") |>
  mutate(region = factor(region, levels = names(region_defs)))

region_cols <- c(
  "MRIP 1" = "violet",  # violet
  "MRIP 2" = "#31A354",  # green
  "MRIP 3" = "red",  # orange
  "MRIP 4" = "#FDD049",  # yellow
  "MRIP 5" = "grey",  # purple-gray
  "MRIP 6" = "#43A2CA",  # blue
  "MRIP 7" = "#F28E2B",  # orange-gold
  "MRIP 8" = "#4b0082"   #  (Keys)
)
pal <- colorFactor(region_cols, levels = names(region_cols), na.color = "#FFFFFF")



m <- leaflet(flc, options = leafletOptions(zoomSnap = 0.25, zoomDelta = 0.25, zoomControl = FALSE)) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  addPolygons(
    color = "black", weight = 1,
    fillColor = ~pal(region), fillOpacity = 0.55,
    label = ~paste0(NAME, ifelse(is.na(region), "", paste0("  (", region, ")"))),
    highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE)
  ) %>%
  addScaleBar(position = "topright") %>%
  addControl(html = "<div style='font-size:20px;font-weight:bold;'>↑<br>N</div>",
             position = "topright") %>%
  addLegend(
    position = "bottomright",
    colors = unname(region_cols),
    labels = names(region_cols),
    opacity = 0.8,
    title = "MRIP Region"
  )

mapview::mapshot(
  m,
  file    = "fl_map.png",
  vwidth  = 2000,
  vheight = 2400,
  zoom    = 2
)
