# Common cartography in a single view, with augmented insular representation.
# Source geometries and departmental counts are preserved.
cargar_cartografia_departamental <- function() {
  dir <- file.path(SIVIGILA_PATH, "datos_externos", "cartografia")
  dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  ruta <- file.path(dir, "dane_mgn2025_departamentos_islas_detalle.rds")
  if (file.exists(ruta)) return(readRDS(ruta))
  original <- file.path(dir, "dane_mgn2025_departamentos.geojson")
  if (!file.exists(original)) {
    url <- paste0("https://geoportal.dane.gov.co/mparcgis/rest/services/",
                  "MGN2025/Serv_CapasMGN_2025/FeatureServer/319/query?",
                  "where=1%3D1&outFields=DPTO_CCDGO%2CDPTO_CNMBRE&returnGeometry=true&outSR=4326&f=geojson")
    utils::download.file(url, original, mode="wb", method="libcurl")
  }
  raw <- sf::st_make_valid(sf::st_read(original, quiet=TRUE))
  raw$codigo_dpto <- sprintf("%02d", as.integer(raw$DPTO_CCDGO))
  raw$nombre_dpto <- raw$DPTO_CNMBRE
  raw <- raw[, c("codigo_dpto", "nombre_dpto")]
  historical <- file.path(dir, "dane_mgn2025_departamentos_simplificado.rds")
  if (file.exists(historical)) {
    mapa <- readRDS(historical)
  } else {
    mapa <- sf::st_transform(sf::st_simplify(sf::st_transform(raw, 9377),
                                             dTolerance=1500, preserveTopology=TRUE), 4326)
  }
  islands <- sf::st_transform(sf::st_simplify(sf::st_transform(raw[raw$codigo_dpto=="88", ], 32617),
                                              dTolerance=20, preserveTopology=TRUE), 4326)
  sf::st_geometry(mapa)[mapa$codigo_dpto=="88"] <- sf::st_geometry(islands)
  stopifnot(nrow(mapa)==33, anyDuplicated(mapa$codigo_dpto)==0,
            all(sf::st_is_valid(mapa)), !any(sf::st_is_empty(mapa)))
  saveRDS(mapa, ruta, compress="xz")
  mapa
}

# Geodesic distance between the ends of a horizontal bar in a metric CRS.
# sf uses the ellipsoid when disabling s2; the previous state is always restored.
distancia_barra_m <- function(x, y, longitud, crs) {
  points <- sf::st_sfc(sf::st_point(c(x,y)), sf::st_point(c(x+longitud,y)), crs=crs)
  geo <- sf::st_transform(points,4326)
  prev <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(prev)))
  as.numeric(sf::st_distance(geo[1],geo[2]))
}

adornar_mapa <- function(g, limites, crs, km, pequeno=FALSE) {
  x0 <- limites[1]; x1 <- limites[2]; y0 <- limites[3]; y1 <- limites[4]
  w <- x1-x0; h <- y1-y0
  # Localized scale calibrated by ellipsoidal distance, not by degrees.
  sx <- x0+.07*w; sy <- y0+.055*h
  long_val <- uniroot(function(l) distancia_barra_m(sx,sy,l,crs)-km*1000,
                      c(km*500, km*1500), tol=.001)$root
  mid_val <- uniroot(function(l) distancia_barra_m(sx,sy,l,crs)-km*500,
                     c(0,long_val), tol=.001)$root
  bh <- h*.012
  g <- g + ggplot2::annotate("rect", xmin=sx, xmax=sx+mid_val, ymin=sy, ymax=sy+bh,
                             fill="#17202A", colour="#17202A", linewidth=.25) +
    ggplot2::annotate("rect", xmin=sx+mid_val, xmax=sx+long_val, ymin=sy, ymax=sy+bh,
                      fill="white", colour="#17202A", linewidth=.25) +
    ggplot2::annotate("text", x=c(sx,sx+mid_val,sx+long_val), y=sy-bh*1.3,
                      label=c("0",scales::number(km/2, accuracy=if(km %% 2) .1 else 1, decimal.mark=","), paste(km,"km")),
                      size=if(pequeno) 2.8 else 3, vjust=1, colour="#17202A")
  # Four-point compass rose oriented to geographic north at its anchor point.
  cx <- x1-.13*w; cy <- y1-.16*h; r <- min(w,h)*.062
  anchor <- sf::st_sfc(sf::st_point(c(cx,cy)), crs=crs)
  ll <- sf::st_coordinates(sf::st_transform(anchor,4326))[1,]
  north <- sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_point(ll+c(0,.01)), crs=4326), crs))[1,]
  theta <- atan2(north[2]-cy,north[1]-cx)
  for (i in 0:3) {
    a <- theta-i*pi/2
    point_val <- c(cx+r*cos(a),cy+r*sin(a))
    left_val <- c(cx+.23*r*cos(a+pi/2),cy+.23*r*sin(a+pi/2))
    right_val <- c(cx+.23*r*cos(a-pi/2),cy+.23*r*sin(a-pi/2))
    g <- g + ggplot2::annotate("polygon", x=c(cx,point_val[1],left_val[1]), y=c(cy,point_val[2],left_val[2]),
                               fill="#17202A", colour="#17202A", linewidth=.2) +
      ggplot2::annotate("polygon", x=c(cx,point_val[1],right_val[1]), y=c(cy,point_val[2],right_val[2]),
                        fill="white", colour="#17202A", linewidth=.2)
  }
  angles <- theta-(0:3)*pi/2
  g <- g + ggplot2::annotate("text", x=cx+1.5*r*cos(angles), y=cy+1.5*r*sin(angles),
                             label=c("N","E","S","W"), fontface="bold", size=if(pequeno) 2.8 else 3)
  attr(g,"control_cartografico") <- list(crs=crs, limites=limites,
                                         escala=c(x=sx,y=sy,longitud=long_val,mitad=mid_val,metros=km*1000),
                                         rosa=c(x=cx,y=cy,angulo_norte=theta), ampliacion=pequeno)
  g
}

dibujar_mapa_departamental <- function(mapa, variable, titulo, subtitulo, caption, etiqueta, n_total) {
  stopifnot(sum(mapa$codigo_dpto=="88")==1L, nrow(mapa)==33L)
  color_limits <- c(0,max(mapa[[variable]], na.rm=TRUE))
  if (color_limits[2]==0) color_limits[2] <- 1
  projected <- sf::st_transform(mapa,9377)
  b <- sf::st_bbox(projected)
  w <- unname(b["xmax"]-b["xmin"]); h <- unname(b["ymax"]-b["ymin"])
  limites <- unname(c(b["xmin"]-.10*w,b["xmax"]+.075*w,b["ymin"]-.10*h,b["ymax"]+.07*h))
  width <- limites[2]-limites[1]; height <- limites[4]-limites[3]
  continent <- mapa[mapa$codigo_dpto!="88",]
  insular <- projected[projected$codigo_dpto=="88",]
  # Integrated and uniformly augmented insular representation. Source geometries
  # and data are preserved; only this drawing layer is scaled.
  insular_factor <- 4
  bb <- sf::st_bbox(insular)
  origin <- unname(c(bb["xmin"],bb["ymin"]))
  destination <- c(limites[1]+.095*width, limites[3]+.81*height)
  orig_geom <- sf::st_geometry(insular)
  pieces <- suppressWarnings(sf::st_cast(orig_geom,"POLYGON"))
  latitude <- sf::st_coordinates(sf::st_transform(sf::st_centroid(pieces),4326))[,2]
  is_providencia <- latitude > 13
  bbox_san <- sf::st_bbox(pieces[!is_providencia])
  bbox_pro <- sf::st_bbox(pieces[is_providencia])
  origin <- unname(c(bbox_san["xmin"],bbox_san["ymin"]))
  origin_pro <- unname(c(bbox_pro["xmin"],bbox_pro["ymin"]))
  # Bring the two groups closer without deforming each island or its orientation. Providencia
  # and Santa Catalina move together; the separation no longer expresses distance.
  destination_pro <- destination+c(45000,unname(bbox_san["ymax"]-bbox_san["ymin"])*insular_factor+40000)
  displacements <- t(vapply(seq_along(pieces),function(i) {
    if(is_providencia[i]) destination_pro-insular_factor*origin_pro else destination-insular_factor*origin
  },numeric(2)))
  polygons <- lapply(seq_along(pieces),function(i) lapply(pieces[[i]],function(ring)
    sweep(ring*insular_factor,2,displacements[i,],"+")))
  transformed <- sf::st_sfc(sf::st_multipolygon(polygons),crs=9377)
  sf::st_geometry(insular) <- transformed
  g <- ggplot2::ggplot(mapa) +
    ggplot2::geom_sf(data=continent,ggplot2::aes(fill=.data[[variable]]),
                     colour="white",linewidth=.20) +
    ggplot2::scale_fill_viridis_c(option="C",trans="sqrt",limits=color_limits,
                                  breaks=scales::pretty_breaks(5),
                                  labels=scales::label_number(big.mark=".",decimal.mark=","),name=etiqueta) +
    ggplot2::labs(title=titulo,subtitle=subtitulo,
                  caption=paste0(caption,
                                 "\nIslands at 4x linear scale; separation adjusted for readability. The island bar does not measure distance between islands.")) +
    ggplot2::theme_void(base_size=11) +
    ggplot2::theme(
      panel.background=ggplot2::element_rect(fill="white",colour=NA),
      plot.background=ggplot2::element_rect(fill="white",colour=NA),
      panel.border=ggplot2::element_rect(fill=NA,colour="#94A3B8",linewidth=.55),
      plot.title=ggplot2::element_text(face="bold",size=14,colour="#17202A",margin=ggplot2::margin(b=5)),
      plot.subtitle=ggplot2::element_text(size=10,colour="#475569",margin=ggplot2::margin(b=10)),
      plot.caption=ggplot2::element_text(size=8,colour="#64748B",hjust=0,lineheight=1.15,margin=ggplot2::margin(t=8)),
      plot.margin=ggplot2::margin(14,14,12,14),
      legend.position="inside",legend.position.inside=c(.96,.06),
      legend.justification.inside=c(1,0),
      legend.background=ggplot2::element_rect(fill="white",colour=NA),
      legend.title=ggplot2::element_text(face="bold",size=11),
      legend.text=ggplot2::element_text(size=10),
      legend.key.height=grid::unit(3.4,"cm"),legend.key.width=grid::unit(.45,"cm")) +
    ggplot2::guides(fill=ggplot2::guide_colourbar(theme=ggplot2::theme(
      legend.key.height=grid::unit(3.4,"cm"),legend.key.width=grid::unit(.45,"cm"))))
  g <- adornar_mapa(g,limites,9377,400)
  control <- attr(g,"control_cartografico")
  g <- g + ggplot2::geom_sf(data=insular,ggplot2::aes(fill=.data[[variable]]),
                            colour="#334155",linewidth=.22,show.legend=FALSE)
  # A compact box gathers the two insular groups and their scale.
  bi <- sf::st_bbox(insular)
  frame <- unname(c(bi["xmin"]-.025*width,bi["xmax"]+.18*width,
                    bi["ymin"]-.065*height,bi["ymax"]+.025*height))
  g <- g + ggplot2::annotate("rect",xmin=frame[1],xmax=frame[2],ymin=frame[3],ymax=frame[4],
                             fill=NA,colour="#94A3B8",linewidth=.5) +
    ggplot2::annotate("text",x=bi["xmax"]+.02*width,y=destination[2]+.011*height,
                      label="San Andrés",hjust=0,size=2.9,colour="#475569") +
    ggplot2::annotate("text",x=bi["xmax"]+.02*width,y=destination_pro[2]+.008*height,
                      label="Providencia and\nSanta Catalina",hjust=0,size=2.9,lineheight=1.1,colour="#475569")
  # 50 km scale for the insular layer: calibrated in its original space
  # and applying exactly the same affine transformation as to the islands.
  ix <- origin[1]; iy <- origin[2]-16000
  long_val <- uniroot(function(l) distancia_barra_m(ix,iy,l,9377)-50000,c(25000,75000),tol=.001)$root
  mid_val <- uniroot(function(l) distancia_barra_m(ix,iy,l,9377)-25000,c(0,long_val),tol=.001)$root
  sx <- destination[1]; sy <- (iy-origin[2])*insular_factor+destination[2]
  bh <- height*.005
  g <- g + ggplot2::annotate("rect",xmin=sx,xmax=sx+mid_val*insular_factor,
                             ymin=sy,ymax=sy+bh,fill="#17202A",colour="#17202A",linewidth=.2) +
    ggplot2::annotate("rect",xmin=sx+mid_val*insular_factor,xmax=sx+long_val*insular_factor,
                      ymin=sy,ymax=sy+bh,fill="white",colour="#17202A",linewidth=.2) +
    ggplot2::annotate("text",x=c(sx,sx+mid_val*insular_factor,sx+long_val*insular_factor),
                      y=sy-bh*2,label=c("0","25","50 km"),size=2.9,vjust=1,colour="#334155")
  # Apply projection after all sf layers, which can incorporate
  # their default coordinate when added to the plot.
  g <- g + ggplot2::coord_sf(crs=sf::st_crs(9377),default_crs=sf::st_crs(9377),
                             xlim=limites[1:2],ylim=limites[3:4],datum=NA,expand=FALSE)
  attr(g,"control_cartografico") <- c(control,list(
    factor_insular=insular_factor,origen_insular=origin,destino_insular=destination,
    desplazamientos_insulares=displacements,marco_insular=frame,
    escala_insular=c(x=ix,y=iy,longitud=long_val,mitad=mid_val,metros=50000),
    n_total=n_total,capas_insulares=c(geometria=14,barra1=18,barra2=19,etiquetas=20)))
  g
}