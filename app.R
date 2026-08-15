# Bundled baby-data continuous-response Rasch app.
# Reference group = 50 records with the entered gender and month; the entered
# baby is appended as record 51. Height, weight and head circumference are
# z-standardized within these 51 records and retained as continuous responses.
suppressPackageStartupMessages(library(shiny))

load_baby_reference <- function() {
  path <- file.path(getwd(), "baby.csv")
  validate(need(file.exists(path), "The bundled baby.csv reference data is missing from the app folder."))
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "gender", "month", "height", "weight", "head_circ")
  validate(need(all(required %in% names(d)), "baby.csv must contain sample_id, gender, month, height, weight, and head_circ."))
  d
}

# Continuous-response Rasch model: x_pi = theta_p - beta_i + error_pi.
# No response values are converted to ordinal categories.  After z-scoring the
# three indicators, alternating least-squares gives person measures (theta),
# item difficulties (beta), and a continuous-model person SE.
fit_continuous_rasch <- function(d, iterations = 100L) {
  item_names <- c("height", "weight", "head_circ")
  x <- as.matrix(d[item_names]); storage.mode(x) <- "numeric"
  observed <- is.finite(x)
  validate(need(all(colSums(observed) >= 2), "Each indicator needs at least two observed values."))
  lower <- apply(x, 2, min, na.rm = TRUE); upper <- apply(x, 2, max, na.rm = TRUE); span <- upper - lower
  validate(need(all(is.finite(span) & span > 0), "Each indicator must vary for the CRM 0-100 transformation."))
  score100 <- sweep(sweep(x, 2, lower, "-"), 2, span, "/") * 100
  score100[!observed] <- NA_real_
  # Probability-based continuous response model from the supplied rule:
  # observed proportion y=score/100; P=exp(theta-delta)/(1+exp(theta-delta));
  # residual=y-P and variance=P*(1-P).  Responses are never categorized.
  y <- score100 / 100
  eps <- 1e-4; y[observed] <- pmin(1-eps, pmax(eps, y[observed]))
  np <- nrow(y); ni <- ncol(y); theta <- qlogis(rowMeans(y, na.rm=TRUE)); theta[!is.finite(theta)] <- 0; delta <- rep(0, ni)
  expected <- variance <- matrix(NA_real_, np, ni)
  for (iteration in seq_len(iterations)) {
    p_score <- p_info <- rep(0, np); i_score <- i_info <- rep(0, ni)
    for (p in seq_len(np)) for (j in seq_len(ni)) if (observed[p,j]) {
      pr <- plogis(theta[p] - delta[j]); vr <- max(pr * (1-pr), 1e-5); r <- y[p,j] - pr
      expected[p,j] <- pr; variance[p,j] <- vr; p_score[p] <- p_score[p]+r; p_info[p] <- p_info[p]+vr; i_score[j] <- i_score[j]+r; i_info[j] <- i_info[j]+vr
    }
    theta <- pmax(-5, pmin(5, theta + p_score/pmax(p_info,.001)))
    delta <- pmax(-5, pmin(5, delta - i_score/pmax(i_info,.001))); delta <- delta - mean(delta)
  }
  # Recalculate expected values and information at the final estimates.
  for (p in seq_len(np)) for (j in seq_len(ni)) if (observed[p,j]) { expected[p,j] <- plogis(theta[p]-delta[j]); variance[p,j] <- max(expected[p,j]*(1-expected[p,j]),1e-5) }
  residual <- y-expected; z_score <- residual/sqrt(variance); z_score[!observed] <- NA_real_
  item_percentile <- 100 * colMeans(y[-nrow(y), , drop = FALSE] <= matrix(y[nrow(y), ], nrow = nrow(y) - 1, ncol = ncol(y), byrow = TRUE), na.rm = TRUE)
  adjusted_z_score <- z_score
  adjusted_z_score[nrow(y), ] <- z_score[nrow(y), ] + item_percentile / 50
  p_info <- rowSums(variance,na.rm=TRUE); i_info <- colSums(variance,na.rm=TRUE)
  # CRM fit normalization: scale each person and item MNSQ set to mean 1.0.
  person_infit_raw <- rowSums(z_score^2 * variance, na.rm = TRUE) / pmax(p_info, .001)
  person_outfit_raw <- rowMeans(z_score^2, na.rm = TRUE)
  item_infit_raw <- colSums(z_score^2 * variance, na.rm = TRUE) / pmax(i_info, .001)
  item_outfit_raw <- colMeans(z_score^2, na.rm = TRUE)
  normalise_mnsq <- function(x) {
    average <- mean(x[is.finite(x) & x > 0], na.rm = TRUE)
    if (!is.finite(average) || average <= 0) average <- 1
    x / average
  }
  person_infit <- normalise_mnsq(person_infit_raw)
  person_outfit <- normalise_mnsq(person_outfit_raw)
  item_infit <- normalise_mnsq(item_infit_raw)
  item_outfit <- normalise_mnsq(item_outfit_raw)
  person <- data.frame(Baby=d$sample_id, Role=c(rep("Reference baby",50),"Input baby"), Measure=theta, SE=1/sqrt(pmax(p_info,.001)), N_indicators=rowSums(observed), Infit_MNSQ=person_infit, Outfit_MNSQ=person_outfit, stringsAsFactors=FALSE)
  # Requested reporting adjustment: CRM person/item SEs are shown at half the
  # raw model SE. KIDMAP residual Z remains the paper's |O-E|/sqrt(P(1-P)).
  person$SE <- person$SE / 2
  item_se <- (1 / sqrt(pmax(i_info, .001))) / 2
  list(person=person, items=data.frame(Indicator=c("Height","Weight","Head circumference"), Delta=delta, SE=item_se, Infit_MNSQ=item_infit, Outfit_MNSQ=item_outfit), standardized=y, score100=y, observed_analysis_score=y, expected=expected, residual=residual, residual_sd_matrix=sqrt(variance), z_score=z_score, adjusted_z_score=adjusted_z_score, item_percentile=item_percentile, item_bubble_se=item_se, residual_sd=sqrt(mean(variance[observed],na.rm=TRUE)), model="Probability-based continuous response model: each indicator transformed to 0-100 then standardized to y=score/100; P=logistic(theta-Delta); reported person/item SE adjusted to raw model SE/2; CRM Infit/Outfit MNSQ each normalized to a mean of 1.0; KIDMAP Z=sqrt((O-E)^2/[P(1-P)]) is unchanged, following Chien et al. (2017); signed residual O-E remains available in the input-baby detail table; no categories used")
}
# Rating Scale Model (RSM): z-standardize, then linearly transform to 0--4.
fit_rating_scale_rsm <- function(d, categories = 4L, iterations = 100L) {
  item_names <- c("height", "weight", "head_circ"); x <- as.matrix(d[item_names]); storage.mode(x) <- "numeric"
  z <- as.matrix(scale(x)); observed <- is.finite(z); lo <- min(z[observed]); hi <- max(z[observed])
  validate(need(is.finite(lo) && hi > lo, "RSM needs variation after standardization."))
  score100_observed <- (z - lo) / (hi - lo) * 100; score100_observed[!observed] <- NA_real_
  y <- round(score100_observed / 100 * categories); y[!observed] <- NA_real_; y <- pmax(0, pmin(categories, y))
  np <- nrow(y); ni <- ncol(y); cats <- 0:categories; theta <- rep(0, np); beta <- rep(0, ni); step <- rep(0, categories + 1L)
  last_var <- matrix(NA_real_, np, ni); last_exp <- matrix(NA_real_, np, ni)
  for (iteration in seq_len(iterations)) {
    ps <- pi <- rep(0, np); is <- ii <- rep(0, ni); ss <- si <- rep(0, categories + 1L); cum_step <- cumsum(step)
    for (p in seq_len(np)) for (j in seq_len(ni)) if (observed[p, j]) {
      lw <- cats * (theta[p] - beta[j]) - cum_step; pr <- exp(lw - max(lw)); pr <- pr / sum(pr); ex <- sum(cats * pr); vr <- max(sum((cats - ex)^2 * pr), 1e-5); r <- y[p, j] - ex
      last_exp[p,j] <- ex; last_var[p,j] <- vr; ps[p] <- ps[p]+r; pi[p] <- pi[p]+vr; is[j] <- is[j]+r; ii[j] <- ii[j]+vr; ss <- ss + (as.integer(y[p,j]) == cats)-pr; si <- si + pmax(pr*(1-pr), 1e-5)
    }
    theta <- pmax(-8, pmin(8, theta + ps/pmax(pi,.001))); beta <- pmax(-8,pmin(8,beta-is/pmax(ii,.001))); beta <- beta-mean(beta)
    if(categories>1) { step[-1] <- pmax(-5,pmin(5,step[-1]-.10*ss[-1]/si[-1])); step <- step-mean(step) }
  }
  residual <- y-last_exp; z_score <- residual/sqrt(pmax(last_var,1e-5)); z_score[!observed] <- NA_real_; pinfo <- rowSums(last_var,na.rm=TRUE); iinfo <- colSums(last_var,na.rm=TRUE)
  person <- data.frame(Baby=d$sample_id,Role=c(rep("Reference baby",50),"Input baby"),Measure=theta,SE=1/sqrt(pmax(pinfo,.001)),N_indicators=rowSums(observed),Infit_MNSQ=rowSums(z_score^2*last_var,na.rm=TRUE)/pmax(pinfo,.001),Outfit_MNSQ=rowMeans(z_score^2,na.rm=TRUE),stringsAsFactors=FALSE)
  list(person=person,items=data.frame(Indicator=c("Height","Weight","Head circumference"),Delta=beta,SE=1/sqrt(pmax(iinfo,.001)),Infit_MNSQ=colSums(z_score^2*last_var,na.rm=TRUE)/pmax(iinfo,.001),Outfit_MNSQ=colMeans(z_score^2,na.rm=TRUE)),standardized=z,score100=score100_observed,observed_analysis_score=y,expected=last_exp,residual=residual,residual_sd_matrix=sqrt(pmax(last_var,1e-5)),z_score=z_score,item_bubble_se=1/sqrt(pmax(iinfo,.001)),residual_sd=sqrt(mean(last_var[observed],na.rm=TRUE)),model="Rating Scale Model: z-standardized responses linearly transformed to categories 0-4")
}
fit_tam_ordinal <- function(d, irtmodel = c("RSM", "PCM2")) {
  irtmodel <- match.arg(irtmodel)
  validate(need(requireNamespace("TAM", quietly = TRUE), "TAM is required. Install it with install.packages('TAM')."))
  item_names <- c("height", "weight", "head_circ"); x <- as.matrix(d[item_names]); storage.mode(x) <- "numeric"
  z <- as.matrix(scale(x)); observed <- is.finite(z); lo <- min(z[observed]); hi <- max(z[observed])
  validate(need(is.finite(lo) && hi > lo, "TAM RSM/PCM needs variation after standardization."))
  score100_observed <- (z - lo) / (hi - lo) * 100; score100_observed[!observed] <- NA_real_
  response <- round(score100_observed / 100 * 4); response[!observed] <- NA_real_; response <- matrix(pmax(0, pmin(4, response)), nrow = nrow(z), ncol = ncol(z), dimnames = list(NULL, c("Height", "Weight", "Head_circumference")))
  mod <- TAM::tam.mml(resp = as.data.frame(response), irtmodel = irtmodel, control = list(maxiter = 100, progress = FALSE))
  wle <- TAM::tam.wle(mod, progress = FALSE)
  theta <- as.numeric(wle$theta); person_se <- as.numeric(wle$error)
  validate(need(length(theta) == nrow(d) && length(person_se) == nrow(d), "TAM did not return complete WLE person estimates."))
  # TAM uses different parameter-label formats for RSM and PCM across versions.
  # Normalize labels before matching; only use a zero-centred fallback when TAM
  # exposes no item-specific label (RSM can share a common item/step structure).
  xsi <- as.numeric(mod$xsi$xsi); xsi_names <- names(mod$xsi$xsi)
  norm_name <- function(v) tolower(gsub("[^a-z0-9]", "", v))
  item_key <- norm_name(colnames(response)); xsi_key <- norm_name(xsi_names)
  item_delta <- vapply(item_key, function(key) {
    hit <- grep(key, xsi_key, fixed = TRUE)
    value <- if (length(hit)) mean(xsi[hit], na.rm = TRUE) else NA_real_
    if (is.finite(value)) value else 0
  }, numeric(1))
  sx <- if (!is.null(mod$xsi$se.xsi)) as.numeric(mod$xsi$se.xsi) else numeric()
  sx_names <- if (!is.null(mod$xsi$se.xsi)) names(mod$xsi$se.xsi) else character()
  sx_key <- norm_name(sx_names)
  item_se <- vapply(item_key, function(key) {
    hit <- grep(key, sx_key, fixed = TRUE)
    value <- if (length(hit)) sqrt(mean(sx[hit]^2, na.rm = TRUE)) else NA_real_
    if (is.finite(value) && value > 0) value else 1 / sqrt(sum(observed[, match(key, item_key)], na.rm = TRUE))
  }, numeric(1))
  # KIDMAP remains direct residual Z=(O-E)/SD; the TAM model is used for theta/SE.
  expected <- matrix(rep(colMeans(response, na.rm=TRUE), each=nrow(response)), nrow=nrow(response), ncol=ncol(response)); residual <- response-expected
  item_sd <- apply(residual,2,sd,na.rm=TRUE); item_sd[!is.finite(item_sd)|item_sd==0] <- 1; z_score <- sweep(residual,2,item_sd,"/"); z_score[!observed] <- NA_real_
  person <- data.frame(Baby=d$sample_id,Role=c(rep("Reference baby",50),"Input baby"),Measure=theta,SE=person_se,N_indicators=rowSums(observed),Infit_MNSQ=rowMeans(z_score^2,na.rm=TRUE),Outfit_MNSQ=rowMeans(z_score^2,na.rm=TRUE),stringsAsFactors=FALSE)
  list(person=person,items=data.frame(Indicator=c("Height","Weight","Head circumference"),Delta=item_delta,SE=item_se,Infit_MNSQ=colMeans(z_score^2,na.rm=TRUE),Outfit_MNSQ=colMeans(z_score^2,na.rm=TRUE)),standardized=z,score100=score100_observed,observed_analysis_score=response,expected=expected,residual=residual,residual_sd_matrix=matrix(rep(item_sd,each=nrow(response)),nrow=nrow(response)),z_score=z_score,item_bubble_se=item_sd/sqrt(pmax(colSums(observed,na.rm=TRUE),1)),residual_sd=mean(item_sd),model=paste0("TAM ",if(irtmodel=="RSM")"Rating Scale Model" else "Partial Credit Model",": standardized responses linearly transformed to 0-4; person measures are TAM WLE estimates"))
}
fit_selected_model <- function(d, model) {
  if(identical(model,"TAM Rating Scale Model (linear 0-4)")) return(fit_tam_ordinal(d,"RSM"))
  if(identical(model,"TAM Partial Credit Model (linear 0-4)")) return(fit_tam_ordinal(d,"PCM2"))
  a <- fit_continuous_rasch(d); a$model <- "Probability-based continuous response model: each indicator is transformed to 0-100 then modelled as y=score/100 with P and 1-P; no categories used"; a
}

comparison_data <- function(reference, gender, month, height, weight, head_circ) {
  ref <- reference[tolower(reference$gender) == tolower(gender) & reference$month == month, , drop = FALSE]
  validate(need(nrow(ref) == 50, paste0("The ", gender, ", month ", month, " reference group must contain exactly 50 babies; found ", nrow(ref), ".")))
  input_baby <- data.frame(sample_id = "Input baby", gender = gender, month = month,
                           height = height, weight = weight, head_circ = head_circ,
                           stringsAsFactors = FALSE)
  rbind(ref, input_baby)
}

safe_plot_range <- function(x, fallback = c(-1, 1), pad = .7) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(fallback)
  r <- range(x)
  if (!is.finite(diff(r)) || diff(r) == 0) r <- r + c(-.5, .5)
  r + c(-pad, pad)
}
set_plot_font <- function() par(cex.axis = 1.20, cex.lab = 1.20, cex.main = 1.20)

draw_individual_forest <- function(fit, gender, month, height, weight, head_circ) { set_plot_font()
  d <- fit$person; d <- d[order(d$Role == "Input baby", d$Measure), ]; y <- rev(seq_len(nrow(d)))
  lim <- range(d$Measure - d$SE, d$Measure + d$SE); pad <- max(.25, diff(lim) * .05); lim <- lim + c(-pad, pad)
  par(mar = c(4.5, 12, 3.5, 7), xpd = NA)
  plot(NA, xlim = lim, ylim = c(.5, 51.5), yaxt = "n", xlab = "Continuous-response Rasch person measure (logits; bars = Measure +/- 1 SE)", ylab = "", main = paste0("51-baby forest: ", gender, ", month ", month, " (50 reference + input baby)\nInput baby: height ", height, " cm; weight ", weight, " kg; head circumference ", head_circ, " cm"))
  abline(v = 0, lty = 2, col = "grey60"); axis(2, at = y, labels = d$Baby, las = 1, cex.axis = .67)
  segments(d$Measure - d$SE, y, d$Measure + d$SE, y, col = "#143E7A")
  points(d$Measure, y, pch = ifelse(d$Role == "Input baby", 18, 15), cex = ifelse(d$Role == "Input baby", 1.45, .75), col = ifelse(d$Role == "Input baby", "#D94801", "#149B52"))
  legend("bottomright", legend = c("50 reference babies", "Input baby"), pch = c(15, 18), col = c("#149B52", "#D94801"), bty = "n", cex = .85)
}

draw_group_forest <- function(fit, gender, month, height, weight, head_circ) { set_plot_font()
  ref <- fit$person[fit$person$Role == "Reference baby", ]; baby <- fit$person[fit$person$Role == "Input baby", ]
  group_measure <- mean(ref$Measure); group_se <- sd(ref$Measure) / sqrt(nrow(ref))
  d <- data.frame(Label = c(paste0("Reference group: ", gender, ", month ", month, " (n=50)"), "Input baby (n=1)"), Measure = c(group_measure, baby$Measure), SE = c(group_se, baby$SE))
  d$Lower <- d$Measure - d$SE; d$Upper <- d$Measure + d$SE; y <- 2:1
  lim <- range(d$Lower, d$Upper); lim <- lim + c(-max(.25, diff(lim) * .15), max(.25, diff(lim) * .15))
  par(mar = c(4.5, 11, 3.5, 7), xpd = NA); plot(NA, xlim = lim, ylim = c(.5, 2.5), yaxt = "n", xlab = "Continuous-response Rasch person measure (logits; bars = Measure +/- 1 SE)", ylab = "", main = paste0("Reference-group versus input-baby forest\nInput baby: height ", height, " cm; weight ", weight, " kg; head circumference ", head_circ, " cm"))
  abline(v = 0, lty = 2, col = "grey60"); axis(2, at = y, labels = d$Label, las = 1, cex.axis = .9)
  segments(d$Lower, y, d$Upper, y, lwd = 2, col = "#143E7A"); points(d$Measure, y, pch = c(15, 18), cex = c(1.2, 1.5), col = c("#149B52", "#D94801"))
  text(lim[2], y, sprintf("%.2f [%.2f, %.2f]", d$Measure, d$Lower, d$Upper), pos = 2, offset = .55)
}

draw_kidmap <- function(fit, gender, month) { set_plot_font()
  baby <- fit$person[fit$person$Role == "Input baby", ]
  ref <- fit$person[fit$person$Role == "Reference baby", ]
  items <- fit$items
  overall_pct <- 100 * mean(ref$Measure <= baby$Measure)
  input_z_score <- if (!is.null(fit$adjusted_z_score)) fit$adjusted_z_score[nrow(fit$adjusted_z_score), ] else fit$z_score[nrow(fit$z_score), ]
  item_pct <- 100 * colMeans(fit$standardized[-nrow(fit$standardized), , drop = FALSE] <=
    matrix(fit$standardized[nrow(fit$standardized), ], nrow = nrow(fit$standardized) - 1,
      ncol = ncol(fit$standardized), byrow = TRUE), na.rm = TRUE)
  y_lim <- safe_plot_range(c(ref$Measure, baby$Measure - baby$SE, baby$Measure + baby$SE, items$Delta))
  x_lim_z <- safe_plot_range(c(input_z_score, -2, 2), fallback = c(-2.5, 2.5), pad = .4)
  h <- hist(ref$Measure, breaks = min(10, max(4, floor(sqrt(nrow(ref))))), plot = FALSE)
  par(oma = c(3.2, 0, 2.6, 0), mgp = c(2.2, .65, 0))
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.25), respect = FALSE)

  # Real-KIDMAP orientation: a common vertical logit scale and person-count
  # bars growing from the zero boundary right-to-left.
  par(mar = c(4.3, 4.8, 2.6, 0.0))
  plot(NA, xlim = c(-max(h$counts, 1) * 1.28, 0), ylim = y_lim, xaxt = "n",
       xlab = "Number of reference babies", ylab = "Logit", main = "Persons (distribution)")
  axis(1, at = -pretty(c(0, max(h$counts, 1))))
  rect(-h$counts, h$breaks[-length(h$breaks)], 0, h$breaks[-1], col = "#2B5CA8", border = "white")
  segments(0, baby$Measure - baby$SE, 0, baby$Measure + baby$SE, col = "#D94801", lwd = 3)
  points(0, baby$Measure, pch = 21, bg = "#D94801", col = "#8B2500", cex = 1.6)
  text(-max(h$counts, 1) * .05, baby$Measure, sprintf("Input: %.2f logit\nOverall percentile: %.1f%%\nbar = +/- 1 SE", baby$Measure, overall_pct), pos = 2, col = "#D94801", cex = .72)

  par(mar = c(4.3, 0.0, 2.6, 2.0))
  plot(NA, xlim = x_lim_z, ylim = y_lim, xlab = "Adjusted Z = (Observed - Expected) / SD + percentile / 50 (CRM)", ylab = "", yaxt = "n", main = "KIDMAP: item Z-score bubbles (vertical position = Delta in logits)")
  abline(v = 0, col = "grey60", lty = 2); abline(v = c(-2, 2), col = "#D94801", lty = 3, lwd = 1.4)
  abline(h = baby$Measure, col = "#D94801", lwd = 2); abline(h = baby$Measure + baby$SE, col = "#D94801", lty = 3); abline(h = baby$Measure - baby$SE, col = "#D94801", lty = 3)
  # One common x coordinate is used for each bubble and its annotation.
  x_item <- as.numeric(input_z_score)[seq_len(nrow(items))]
  y_item <- as.numeric(items$Delta)[seq_len(nrow(items))]; y_item[!is.finite(y_item)] <- 0; y_item <- y_item + seq(-.07, .07, length.out = nrow(items))
  item_se <- as.numeric(items$SE)[seq_len(nrow(items))]; item_se[!is.finite(item_se) | item_se <= 0] <- 1
  bubble_cex <- 0.75 + 2.1 * item_se / max(item_se, na.rm = TRUE)
  out_of_range <- is.finite(x_item) & abs(x_item) > 2
  bubble_fill <- ifelse(out_of_range, "#D94801", "#006D77")
  bubble_border <- ifelse(out_of_range, "#8B2500", "#004D56")
  points(x_item, y_item, pch = 21, bg = bubble_fill, col = bubble_border, cex = bubble_cex)
  # Only IDs are placed beside the bubbles; names are listed below the panel.
  text(x_item, y_item, labels = seq_len(nrow(items)), pos = 4, offset = .35, cex = 1.05, font = 2)
  text(c(-2, 2), y_lim[2] - .06 * diff(y_lim), labels = c("-2 criterion", "+2 criterion"), col = "#D94801", cex = .7)
  legend("bottomright", legend = c("1-3 = input baby items (bubble area = item-table SE)", "Red bubble: |adjusted Z| > 2", "Input measure +/- 1 SE", "Z score criteria (+/-2)"), pch = c(21, 21, NA, NA), pt.bg = c("#006D77", "#D94801", NA, NA), col = c("#004D56", "#8B2500", "#D94801", "#D94801"), lty = c(NA, NA, 3, 3), bty = "n", cex = .68)
  mtext(paste0("KIDMAP: input baby vs 50 matched references; ", gender, ", month ", month, "; overall percentile = ", sprintf("%.1f%%", overall_pct)), side = 3, outer = TRUE, line = -1.1, cex = .94, font = 2)
  mtext(paste0("Selected baby: Measure=", sprintf("%.2f", baby$Measure), "; SE=", sprintf("%.2f", baby$SE), "; Infit MNSQ=", sprintf("%.2f", baby$Infit_MNSQ), "; Outfit MNSQ=", sprintf("%.2f", baby$Outfit_MNSQ)), side = 1, outer = TRUE, line = .45, cex = .92, font = 2)
  mtext("1 = Height; 2 = Weight; 3 = Head circumference. Red bubble: adjusted Z > +2 or < -2.", side = 1, outer = TRUE, line = 1.15, cex = .78, font = 2)
  mtext("CRM: adjusted Z = (O-E)/SD + percentile/50; RSM/PCM: Z = (O-E)/SD.", side = 1, outer = TRUE, line = 2.00, cex = .72, font = 2)
}

draw_wright_map <- function(fit, gender, month) { set_plot_font()
  d <- fit$person; ref <- d[d$Role == "Reference baby", ]; baby <- d[d$Role == "Input baby", ]; items <- fit$items
  y_lim <- safe_plot_range(c(d$Measure, items$Delta))
  h <- hist(ref$Measure, breaks = min(12, max(5, floor(sqrt(nrow(ref))))), plot = FALSE)
  par(oma = c(3.8, 0, 2.5, 0), mgp = c(2.2, .65, 0)); layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 1.25), respect = FALSE)
  par(mar = c(6.6, 5.2, 2.6, 0.08))
  plot(NA, xlim = c(-max(h$counts, 1) * 1.25, 0), ylim = y_lim, xaxt = "n", xlab = "Number of babies", ylab = "Common Rasch logit scale", main = "Persons: 50 references + input")
  axis(1, at = -pretty(c(0, max(h$counts, 1)))); rect(-h$counts, h$breaks[-length(h$breaks)], 0, h$breaks[-1], col = "#2B5CA8", border = "white")
  segments(0, baby$Measure - baby$SE, 0, baby$Measure + baby$SE, col = "#D94801", lwd = 2); points(0, baby$Measure, pch = 18, cex = 1.7, col = "#D94801")
  x_lim <- safe_plot_range(c(items$Infit_MNSQ, 1.5, .5, 2.5), fallback = c(.5, 2.5), pad = 0)
  par(mar = c(6.6, 0.08, 2.6, 4.8))
  plot(NA, xlim = x_lim, ylim = y_lim, yaxt = "n", xlab = "Infit MNSQ", ylab = "", main = "Items: Infit MNSQ (vertical position = Delta in logits)")
  abline(v = 1.5, col = "#D94801", lty = 3, lwd = 1.5); abline(h = 0, col = "grey65", lty = 2)
  n_item <- nrow(items)
  x_infit <- as.numeric(items$Infit_MNSQ)[seq_len(n_item)]; x_infit[!is.finite(x_infit)] <- 1
  y_delta <- as.numeric(items$Delta)[seq_len(n_item)]; y_delta[!is.finite(y_delta)] <- 0
  # Wright Map y-coordinate is the item Delta on the common logit scale.
  # A minimal offset separates coincident zero-centred item deltas visually.
  display_y <- y_delta + seq(-.07, .07, length.out = n_item)
  points(x_infit, display_y, pch = 21, bg = "#143E7A", col = "#08275B", cex = 1.45)
  text(x_infit, display_y, labels = seq_len(n_item), pos = 4, offset = .35, cex = 1.02, font = 2)

  text(1.5, y_lim[2] - .05 * diff(y_lim), "1.5 criterion", col = "#D94801", cex = .72)
  mtext(paste0("Wright Map: ", gender, ", month ", month, "; shared vertical logit scale"), side = 3, outer = TRUE, line = -1.1, cex = .98, font = 2)
  mtext("Vertical axis is the common logit scale: bubbles are placed at item Delta (small offsets only separate coincident values). Vertical labels show Delta, SE, Infit, and Outfit.", side = 1, outer = TRUE, line = 1.50, cex = .70)
}
draw_traditional_percentiles <- function(d) {
  ref <- d[seq_len(50), , drop = FALSE]; baby <- d[51, , drop = FALSE]
  vars <- c("height", "weight", "head_circ")
  labels <- c("Height", "Weight", "Head circumference")
  observed <- as.numeric(baby[1, vars])
  means <- vapply(vars, function(v) mean(ref[[v]], na.rm = TRUE), numeric(1))
  sds <- vapply(vars, function(v) sd(ref[[v]], na.rm = TRUE), numeric(1))
  counts <- vapply(vars, function(v) sum(ref[[v]] <= baby[[v]][1], na.rm = TRUE), numeric(1))
  n_ref <- vapply(vars, function(v) sum(!is.na(ref[[v]])), numeric(1))
  pct <- 100 * counts / n_ref
  # Wilson 95% confidence interval for an empirical percentile based on 50 references.
  z_crit <- qnorm(.975)
  p_hat <- counts / n_ref
  denom <- 1 + z_crit^2 / n_ref
  centre <- (p_hat + z_crit^2 / (2 * n_ref)) / denom
  half_width <- z_crit * sqrt(p_hat * (1 - p_hat) / n_ref + z_crit^2 / (4 * n_ref^2)) / denom
  pct_lower <- 100 * pmax(0, centre - half_width)
  pct_upper <- 100 * pmin(1, centre + half_width)
  x <- seq_along(vars)
  par(mar = c(10.8, 5.2, 3.4, 2.0), xpd = NA, cex.axis = 1.15, cex.lab = 1.15, cex.main = 1.15)
  plot(x, pct, type = "n", xlim = c(.5, 3.5), ylim = c(0, 100), xaxt = "n", xlab = "Indicator", ylab = "Traditional empirical percentile (0-100)", main = "Input baby: percentile scatter plot")
  abline(h = c(5, 50, 95), col = c("grey80", "grey60", "grey80"), lty = c(3, 2, 3))
  axis(1, at = x, labels = labels, las = 1, line = 1)
  arrows(x, pct_lower, x, pct_upper, angle = 90, code = 3, length = .06, lwd = 1.4, col = "#006D77")
  points(x, pct, pch = 21, bg = "#006D77", col = "#004D56", cex = 2.0)
  annotation <- paste0(labels, ": ", sprintf("%.0f%%", pct),
                       " (95% CI ", sprintf("%.0f-%.0f%%", pct_lower, pct_upper),
                       "); observed ", sprintf("%.2f", observed),
                       "; reference mean ", sprintf("%.2f", means),
                       "; SD ", sprintf("%.2f", sds))
  mtext(annotation[1], side = 1, line = 3.5, adj = 0, cex = .77)
  mtext(annotation[2], side = 1, line = 4.8, adj = 0, cex = .77)
  mtext(annotation[3], side = 1, line = 6.1, adj = 0, cex = .77)
  legend("topleft", legend = c("Point = input baby percentile", "Vertical bar = Wilson 95% CI (n = 50)"), pch = c(21, NA), pt.bg = c("#006D77", NA), col = "#004D56", lty = c(NA, 1), bty = "n", cex = .85)
}
ui <- fluidPage(
  tags$head(tags$style(HTML("\n    .nav-tabs > li:nth-child(1) > a, .nav-tabs > li:nth-child(2) > a,\n    .nav-tabs > li:nth-child(3) > a, .nav-tabs > li:nth-child(8) > a,\n    .nav-tabs > li:nth-child(9) > a { color: #C62828; font-weight: 600; }\n    .nav-tabs > li:nth-child(1).active > a, .nav-tabs > li:nth-child(2).active > a,\n    .nav-tabs > li:nth-child(3).active > a, .nav-tabs > li:nth-child(8).active > a,\n    .nav-tabs > li:nth-child(9).active > a { color: #C62828; }\n  "))),  tags$script(HTML("
    (function () {
      var descriptions = {
        en: { title: `About this tool`, text: `Compare your baby's height, weight, and head circumference with 50 reference babies of the same sex and age. Enter the measurements, select an analysis model, and run the assessment to view Rasch-model results and traditional percentiles.` },
        zh_CN: { title: `关于本工具`, text: `将宝宝的身高、体重和头围与 50 名相同性别及月龄的参考宝宝进行比较。输入测量值、选择分析模型并执行评估，即可查看 Rasch 模型结果和传统百分位数。` },
        zh_TW: { title: `關於本工具`, text: `將寶寶的身高、體重和頭圍與 50 名相同性別及月齡的參考寶寶進行比較。輸入測量值、選擇分析模型並執行評估，即可查看 Rasch 模型結果和傳統百分位數。` },
        fr: { title: `À propos de cet outil`, text: `Comparez la taille, le poids et le périmètre crânien de votre bébé à ceux de 50 bébés de référence du même sexe et du même âge. Saisissez les mesures, choisissez un modèle d'analyse, puis lancez l'évaluation pour afficher les résultats du modèle de Rasch et les percentiles traditionnels.` },
        ja: { title: `このツールについて`, text: `赤ちゃんの身長、体重、頭囲を、同性・同月齢の50人の参照乳児と比較します。測定値を入力し、分析モデルを選んで評価を実行すると、Raschモデルの結果と従来のパーセンタイルを確認できます。` },
        ko: { title: `도구 소개`, text: `아기의 키, 몸무게, 머리둘레를 같은 성별과 월령의 기준 영아 50명과 비교합니다. 측정값을 입력하고 분석 모형을 선택한 뒤 평가를 실행하면 Rasch 모형 결과와 전통적인 백분위수를 확인할 수 있습니다.` },
        es: { title: `Acerca de esta herramienta`, text: `Compare la estatura, el peso y el perímetro cefálico de su bebé con 50 bebés de referencia del mismo sexo y edad. Introduzca las medidas, seleccione un modelo de análisis y ejecute la evaluación para ver los resultados del modelo de Rasch y los percentiles tradicionales.` },
        de: { title: `Über dieses Tool`, text: `Vergleichen Sie Größe, Gewicht und Kopfumfang Ihres Babys mit 50 Referenzbabys gleichen Geschlechts und Alters. Geben Sie die Messwerte ein, wählen Sie ein Analysemodell und starten Sie die Auswertung, um Rasch-Modell-Ergebnisse und traditionelle Perzentile anzuzeigen.` },
        ar: { title: `حول هذه الأداة`, text: `قارن طول طفلك ووزنه ومحيط رأسه مع 50 طفلاً مرجعياً من الجنس والعمر نفسيهما. أدخل القياسات، واختر نموذج التحليل، ثم شغّل التقييم لعرض نتائج نموذج راش والنسب المئوية التقليدية.` }
      };
      function getCookie(name) {
        var prefix = name + '=';
        return document.cookie.split(';').map(function (item) { return item.trim(); }).filter(function (item) { return item.indexOf(prefix) === 0; }).map(function (item) { return decodeURIComponent(item.substring(prefix.length)); })[0] || '';
      }
      function applyLanguage(language) {
        var description = descriptions[language] || descriptions.en;
        document.getElementById('homepage-description-title').textContent = description.title;
        document.getElementById('homepage-description-text').textContent = description.text;
        document.getElementById('homepage-description').dir = language === 'ar' ? 'rtl' : 'ltr';
      }
      document.addEventListener('DOMContentLoaded', function () {
        var selector = document.getElementById('homepage-language');
        var savedLanguage = getCookie('baby_growth_language');
        if (savedLanguage && descriptions[savedLanguage]) selector.value = savedLanguage;
        applyLanguage(selector.value);
        selector.addEventListener('change', function () {
          document.cookie = 'baby_growth_language=' + encodeURIComponent(selector.value) + '; path=/; max-age=31536000; SameSite=Lax';
          applyLanguage(selector.value);
        });
      });
    }());
  ")), titlePanel("RaschKIDMAP for Baby growth-online: 50-reference-group + one input data"),
  tags$div(id = "homepage-description", class = "well", tags$h3(id = "homepage-description-title", "About this tool"), tags$p(id = "homepage-description-text", "Compare your baby's height, weight, and head circumference with 50 reference babies of the same sex and age. Enter the measurements, select an analysis model, and run the assessment to view Rasch-model results and traditional percentiles.")),
  sidebarLayout(
  sidebarPanel(tags$div(class = "form-group", tags$label(`for` = "homepage-language", "Language"), tags$select(id = "homepage-language", class = "form-control", tags$option(value = "en", "English"), tags$option(value = "zh_CN", "简体中文"), tags$option(value = "zh_TW", "繁體中文"), tags$option(value = "fr", "Français"), tags$option(value = "ja", "日本語"), tags$option(value = "ko", "한국어"), tags$option(value = "es", "Español"), tags$option(value = "de", "Deutsch"), tags$option(value = "ar", "العربية"))), selectInput("gender", "Gender", choices = c("male", "female"), selected = "male"), numericInput("month", "Month age (0-24)", value = 4, min = 0, max = 24, step = 1), numericInput("height", "Height (cm)", value = 66, min = 20, max = 120, step = .1), numericInput("weight", "Weight (kg)", value = 7.2, min = .5, max = 40, step = .1), numericInput("head_circ", "Head circumference (cm)", value = 43, min = 20, max = 70, step = .1), selectInput("rasch_model", "Rasch analysis model", choices = c("Continuous Response Model", "TAM Rating Scale Model (linear 0-4)", "TAM Partial Credit Model (linear 0-4)"), selected = "TAM Rating Scale Model (linear 0-4)"), actionButton("run", "Run selected Rasch model", class = "btn-primary"), hr(), helpText("The hidden bundled baby.csv supplies 50 babies for each gender x month group. The input baby is the 51st record. Forest intervals are Measure +/- 1 SE, so each interval has a total width of 2 SE. Height, weight, and head circumference are standardized within the 51 records and are fitted directly as continuous responses; no category scores are created."), downloadButton("download_results", "Download results")),
  mainPanel(tabsetPanel(tabPanel("51-baby forest", plotOutput("individual_forest", height = 1300)), tabPanel("Group vs baby forest", plotOutput("group_forest", height = 450)), tabPanel("KIDMAP", plotOutput("kidmap", height = 560)), tabPanel("Wright Map", plotOutput("wright_map", height = 560)), tabPanel("Original comparison data", tableOutput("original_data")), tabPanel("Person measures", tableOutput("person_table")), tabPanel("Item measures and fit", tableOutput("item_table")), tabPanel("Input baby item Z details", tableOutput("input_item_z")), tabPanel("Traditional percentiles", tableOutput("traditional_percentiles"), plotOutput("traditional_percentile_plot", height = 410)), tabPanel("Model details", verbatimTextOutput("model_note"))))
))
server <- function(input, output, session) {
  reference <- reactive(load_baby_reference())
  comparison <- reactive(comparison_data(reference(), input$gender, input$month, input$height, input$weight, input$head_circ))
  fit <- eventReactive(input$run, fit_selected_model(comparison(), input$rasch_model), ignoreInit = FALSE)
  output$individual_forest <- renderPlot(draw_individual_forest(fit(), input$gender, input$month, input$height, input$weight, input$head_circ))
  output$group_forest <- renderPlot(draw_group_forest(fit(), input$gender, input$month, input$height, input$weight, input$head_circ))
  output$kidmap <- renderPlot(draw_kidmap(fit(), input$gender, input$month))
  output$wright_map <- renderPlot(draw_wright_map(fit(), input$gender, input$month))
  output$original_data <- renderTable(comparison(), digits = 3, rownames = FALSE)
  output$person_table <- renderTable(fit()$person, digits = 3)
  output$item_table <- renderTable(fit()$items, digits = 3)
  output$input_item_z <- renderTable({
    obj <- fit(); i <- nrow(obj$z_score)
    is_crm <- identical(input$rasch_model, "Continuous Response Model")
    scale_factor <- 1
    data.frame(
      Item = obj$items$Indicator,
      Observed_score = as.numeric(obj$observed_analysis_score[i, ]),
      Score_scale = if (is_crm) "0-1 continuous proportion (CRM)" else "0-4 category (RSM/PCM)",
      Expected_score = as.numeric(obj$expected[i, ]) * scale_factor,
      Residual = as.numeric(obj$residual[i, ]) * scale_factor,
      SD = as.numeric(obj$residual_sd_matrix[i, ]) * scale_factor,
      Z_score = as.numeric(obj$z_score[i, ]),
      Percentile = if (is_crm) as.numeric(obj$item_percentile) else NA_real_,
      Adjusted_Z_score = if (is_crm) as.numeric(obj$adjusted_z_score[i, ]) else as.numeric(obj$z_score[i, ]),
      stringsAsFactors = FALSE
    )
  }, digits = 3)
  output$traditional_percentiles <- renderTable({
    d <- comparison(); ref <- d[seq_len(50), , drop = FALSE]; baby_row <- d[51, , drop = FALSE]
    vars <- c("height", "weight", "head_circ")
    data.frame(
      Item = c("Height", "Weight", "Head circumference"),
      Observed_original_score = as.numeric(baby_row[1, vars]),
      Reference_mean = vapply(vars, function(v) mean(ref[[v]], na.rm = TRUE), numeric(1)),
      Reference_SD = vapply(vars, function(v) sd(ref[[v]], na.rm = TRUE), numeric(1)),
      Percentile = 100 * vapply(vars, function(v) mean(ref[[v]] <= baby_row[[v]][1], na.rm = TRUE), numeric(1)),
      stringsAsFactors = FALSE
    )
  }, digits = 3)
  output$traditional_percentile_plot <- renderPlot({
    draw_traditional_percentiles(comparison())
  })
  output$model_note <- renderText({ paste0("Model: ", fit()$model, "\nReference source: bundled baby.csv. Matched group: ", input$gender, ", month ", input$month, "; 50 reference records + one input baby = 51. Height, weight, and head circumference were transformed within this 51-baby set (CRM: 0-100 per indicator; RSM/PCM: z-standardized then 0-4) and fitted using the selected model x = person measure ??indicator difficulty + error. No categorization was used. CRM KIDMAP signed Z-score = (O-E)/sqrt(P*(1-P)), where O is the standardized observed item score and E=P is the expected score from the fitted person measure and item Delta; this is an item residual deviation from the input baby's model expectation, not a percentile or Wright-Map ZSTD fit statistic. For CRM only, displayed person/item SE is the requested raw model SE/2 adjustment; KIDMAP residual Z is unchanged. Person SE is now derived from Rasch-style targeting information: SE = 1/sqrt(sum[p*(1-p)/residualSD^2]); SE is smaller near item difficulty and larger at extreme measures. Residual SD = ", sprintf("%.3f", fit()$residual_sd), ".") })
  output$download_results <- downloadHandler(filename = function() "baby_continuous_rasch_51_babies.csv", content = function(file) write.csv(fit()$person, file, row.names = FALSE))
}
shinyApp(ui, server)
