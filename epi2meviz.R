library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(plotly)
library(shinydashboard)
library(shinyBS)
library(DT)
library(htmlwidgets)
library(leaflet)
library(tidyr)
library(shinyjs)
library(microbiome)
library(phyloseq)
library(vegan)
library(shinyalert)

### Load Data ###
ui <- dashboardPage(skin = "red",
  dashboardHeader(title = "EPI2MEViz - Reboot"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("file")),
      menuItem("Rarefaction", tabName = "rarefaction", icon = icon("chart-area")),
      menuItem("Relative Abundance", tabName = "relativeabundance", icon = icon("chart-pie")),
      menuItem("Alpha Diversity", tabName = "alphadiversity", icon = icon("chart-bar")),
      menuItem("PCoA", tabName = "pcoa", icon = icon("chart-line"))
    )
  ),

  # Formatting to make it iolani colors
  dashboardBody(tags$style(HTML("
  .box.box-solid.box-primary>.box-header {
    color:#fff;
    background:#000000
  }

  .box.box-solid.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .box.box-primary>.box-header {
    color:#000000;
    background:#fff
  }

  .box.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .skin-red .main-sidebar {
    background-color: #000000;
  }
  ")),

  tabItems(
    tabItem(tabName = "upload",
      fluidRow(
        box(
          title = "EPI2ME Taxa Table Upload", status = "primary", solidHeader = TRUE, width = 12,
          fileInput("data", "Input csv files", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"), multiple = FALSE),
          actionButton("datasubmit", "Submit"),
          actionButton("datahelp", "Help!")
        )
      ),
      fluidRow(
        box(
          title = "Metadata Upload (Optional)", status = "primary", solidHeader = TRUE, width = 12,
          fileInput("meta", "Input csv files", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"), multiple = FALSE),
          actionButton("metadatahelp", "Help!")
        )
      )
    ),
    tabItem(tabName = "rarefaction",
      fluidRow(
        box(
          style = "height: 100%;",
          title = "Rarefaction Curve", status = "primary", solidHeader = TRUE, width = 12,
          plotlyOutput("rarecurveplot",
            height = "100%"
          ),
          actionButton("rarehelp", "Help!")
        )
      )
    ),
    tabItem(tabName = "relativeabundance",
      fluidRow(
        box(
          style = "height: 100%;",
          title = "Relative Abundance", status = "primary", solidHeader = TRUE, width = 12,
          selectInput("taxlevel1", "Select Taxonomic Level:", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), selected = "Genus"),
          plotlyOutput("relabundplot",
            height = "100%"
          ),
          actionButton("relabundhelp", "Help!")
        )
      )
    ),
    tabItem(tabName = "alphadiversity",
      fluidRow(
        box(
          style = "height: 100%;",
          title = "Alpha Diversity", status = "primary", solidHeader = TRUE, width = 12,
          selectInput("grouping", "Select Grouping Variable (if metadata provided):", choices = c("barcode"), selected = "barcode"),
          plotlyOutput("alphadivplot",
            height = "100%"
          ),
          actionButton("alphadivhelp", "Help!")
        )
      )
    ),
    tabItem(tabName = "pcoa",
      fluidRow(
        box(
          style = "height: 100%;",
          title = "Principal Coordinate Analysis", status = "primary", solidHeader = TRUE, width = 12,
          selectInput("taxlevel2", "Select Taxonomic Level:", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), selected = "Genus"),
          selectInput("colorby", "Select Variable to Color By (if metadata provided):", choices = c("barcode"), selected = "barcode"),
          plotlyOutput("pcoaplot",
            height = "100%"
          ),
          actionButton("pcoahelp", "Help!")
        )
      )
    )
  ))
)

server <- function(input, output, session) {
  # Set large upload size limit (server side)
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Load help modules
  source("modules/help.R", local = TRUE)

  dat <- reactiveVal(NULL)

  observeEvent(input$datasubmit, {
    # Require data upload before proceeding, metadata is optional
    req(input$data)

    # Ingest data
    df <- read.csv(input$data$datapath, header = TRUE, stringsAsFactors = FALSE)
    message("Uploaded data: ", nrow(df), " rows, ", ncol(df), " cols")

    # Combine data and metadata if applicable
    if (!is.null(input$meta)) {
      # Process with metadata
      req(input$meta)
      meta_df <- read.csv(input$meta$datapath, header = TRUE, stringsAsFactors = FALSE)
      # Coerce all metdata to character
      meta_df <- meta_df %>%
        mutate(across(everything(), as.character))

      # Convert all metadata columns to lower
      colnames(meta_df) <- tolower(colnames(meta_df))

      message("Uploaded metadata: ", nrow(meta_df), " rows, ", ncol(meta_df), " cols")

      # Extract columns with "barcode" in their names
      otu_dat <- df %>%
        select(starts_with("barcode"))
      otu_dat <- otu_table(otu_dat, taxa_are_rows = TRUE)

      # Extract taxonomic information
      tax_dat <- df %>%
        select(superkingdom, phylum, class, order, family, genus)

      tax_dat <- tax_table(tax_dat)
      taxa_names(otu_dat) <- rownames(tax_dat)

      # Drop metadata rows without matching sample names in otu_dat
      meta_df <- meta_df %>%
        filter(barcode %in% colnames(otu_dat))

      # Ensure metadata rownames match sample names in otu_dat
      rownames(meta_df) <- meta_df$barcode

      # Create phyloseq object
      physeq <- phyloseq(otu_dat, tax_dat, sample_data(meta_df))

      # Name the columns of phyloseq object
      colnames(physeq@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      print(physeq)
      dat(physeq)
    } else {
      # Process without metadata
      message("No metadata uploaded.")

      # Extract columns with "barcode" in their names
      otu_dat <- df %>%
        select(starts_with("barcode"))
      otu_dat <- otu_table(otu_dat, taxa_are_rows = TRUE)

      # Extract taxonomic information
      tax_dat <- df %>%
        select(superkingdom, phylum, class, order, family, genus)

      tax_dat <- tax_table(tax_dat)
      taxa_names(otu_dat) <- rownames(tax_dat)

      # Create metdata placeholder with barcodes as rownames
      # Needed to add an extra dummy column to avoid PCoA color issues
      meta <- data.frame(colnames(otu_dat), colnames(otu_dat))
      rownames(meta) <- colnames(otu_dat)
      names(meta) <- c("barcode", "barcode_extra")

      # Create phyloseq object
      physeq <- phyloseq(otu_dat, tax_dat, sample_data(meta))

      # Name the columns of phyloseq object
      colnames(physeq@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      print(physeq)
      dat(physeq)
    }

    # Show modal to indicate successfull analysis
    shinyalert::shinyalert(
      title = "Data Uploaded Successfully!",
      text = paste("The data has been uploaded and processed. You can now navigate to the other tabs to explore the visualizations."),
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK"
    )
  })

  # Update grouping (and colorby) selectInput choices from metadata columns when data is available
  observeEvent(dat(), {
    req(dat())
    # Grab current phyloseq object and its sample metadata
    phy <- dat()
    md <- as.data.frame(sample_data(phy))
    cols <- colnames(md)
    if (length(cols) == 0) {
      cols <- "barcode"
    }

    # Prefer a sensible default: first metadata column that isn't a barcode helper
    non_barcode <- setdiff(cols, c("barcode", "barcode_extra"))
    sel <- if (length(non_barcode) > 0) non_barcode[1] else cols[1]

    updateSelectInput(session, "grouping", choices = cols, selected = sel)
    # also update the color selector used in other plots
    updateSelectInput(session, "colorby", choices = cols, selected = sel)

    # Convert non-numeric metadata columns to factors so grouping works correctly in plots
    if (ncol(md) > 0) {
      md[] <- lapply(md, function(x) if (!is.numeric(x)) as.factor(x) else x)
      sample_data(phy) <- sample_data(md)
      dat(phy)
    }
  })

  # Rarefaction curve generation
  output$rarecurveplot <- renderPlotly({
    req(dat())
    raref <- dat()

    raref <- prune_taxa(taxa_sums(raref) > 1, raref)
    otu_rare <- raref@otu_table
    otu_rare <- as.data.frame(otu_rare)
    sample_names <- rownames(otu_rare)

    rareplot <- rarecurve(t(otu_rare), step = 50, sample = 5000, tidy = TRUE)

    ggplotly(
      ggplot(rareplot, aes(x = Sample, y = Species, color = Site)) +
        geom_line() +
        theme_classic() +
        labs(x = "Sequencing Depth", y = "Observed Species")
    )
  })

  # Relative abundance plot generation
  output$relabundplot <- renderPlotly({
    req(dat())
    relab <- dat()
    relab <- tax_glom(relab, taxrank = input$taxlevel1)
    relab <- transform_sample_counts(relab, function(x) x / sum(x) * 100)

    ggplotly(
      plot_bar(relab, fill = input$taxlevel1, title = "Relative Abundance Plot") +
        geom_bar(stat = "identity", position = "stack", color = NA) +
        theme_classic() +
        ylab("Relative Abundance") +
        xlab("Samples") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
    )
  })

  # Alpha diversity plot generation
  output$alphadivplot <- renderPlotly({
    req(dat())
    alpha_div <- dat()

    print(alpha_div@sam_data)

    ggplotly(
      plot_richness(alpha_div, x = input$grouping, measures = c("Simpson"), color = input$grouping) +
        geom_boxplot() +
        theme_classic() +
        xlab("Samples") +
        ylab("Diversity Index") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), legend.position = "none")
    )
  })

  # PCoA
  output$pcoaplot <- renderPlotly({
    req(dat())
    pcoa_dat <- dat()
    pcoa_dat <- tax_glom(pcoa_dat, taxrank = input$taxlevel2)

    # Create a distance matrix based on bray-curtis dissimilarity
    ord <- ordinate(pcoa_dat, "MDS", "bray")

    plot_ordination(pcoa_dat, ord, color = input$colorby) +
      geom_point(size = 5) +
      theme_classic()
  })

  #####################
  ### HELP MESSAGES ###
  #####################
  observeEvent(input$datahelp, {
    shinyalert(
      title = "Data Upload Help",
      text = paste0("To upload your data, please ensure you have a CSV file with OTU counts and taxonomic information. The OTU counts should be in columns starting with 'barcode', and taxonomic information should include columns for superkingdom, phylum, class, order, family, and genus. Optionally, you can upload a separate metadata CSV file with sample information."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })

  observeEvent(input$metadatahelp, {
    shinyalert(
      title = "Metadata Upload Help",
      text = paste0("Uploading metadata is optional but recommended for enhanced analysis. Your metadata CSV file should contain a 'barcode' column that matches the sample names in your OTU table. Additional columns can include any relevant sample information you wish to use for grouping or coloring in visualizations."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })

  observeEvent(input$rarehelp, {
    shinyalert(
      title = "Rarefaction Curve Help",
      text = paste0("The rarefaction curve visualizes the relationship between sequencing depth and observed species richness. It helps assess whether the sequencing effort was sufficient to capture the diversity of the microbial community. Each line represents a sample, showing how many species were observed at different levels of sequencing depth."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })

  observeEvent(input$relabundhelp, {
    shinyalert(
      title = "Relative Abundance Plot Help",
      text = paste0("The relative abundance plot displays the proportion of different taxa at the selected taxonomic level across samples. It helps visualize the composition of microbial communities and identify dominant taxa. You can select the taxonomic level (e.g., Genus, Family) to explore how microbial populations vary among samples."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })

  observeEvent(input$alphadivhelp, {
    shinyalert(
      title = "Alpha Diversity Plot Help",
      text = paste0("The alpha diversity plot illustrates the diversity within individual samples using the selected diversity index (e.g., Simpson). It helps compare the richness and evenness of microbial communities across different groups or conditions. You can select a grouping variable from your metadata to see how diversity varies among different sample categories."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })

  observeEvent(input$pcoahelp, {
    shinyalert(
      title = "PCoA Plot Help",
      text = paste0("The Principal Coordinate Analysis (PCoA) plot visualizes the differences in microbial community composition between samples based on a distance matrix (e.g., Bray-Curtis dissimilarity). It helps identify patterns and clusters among samples. You can select a taxonomic level for analysis and a variable from your metadata to color the points, allowing for better interpretation of the relationships between samples."),
      type = "info",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })
}

shinyApp(ui, server)