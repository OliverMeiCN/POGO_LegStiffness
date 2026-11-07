# ============================================================
# POGO
# Automated Leg Stiffness Analysis
# Task: Countermovement Jump (CMJ) Landing
#
# Software name: POGO
# Module: CMJ landing leg stiffness analysis from C3D files
#
# The original event-detection and stiffness-calculation logic
# is retained. Only the software identification and output naming
# have been updated for the POGO project.
# ============================================================

if (!requireNamespace("c3dr", quietly = TRUE)) {
  install.packages("c3dr", dependencies = TRUE)
}

library(c3dr)

cat("\n========== POGO | CMJ Leg Stiffness Analysis ==========\n")
cat("Software: POGO\nTask: CMJ landing\n")
cat("当前 c3dr 版本：\n")
print(packageVersion("c3dr"))


# ============================================================
# 1. 设置要批量处理的文件夹
# ============================================================

input_folder <- "E:/stiffness/leg stiffness/data/leg stiffness_data/jump_male/addition/add"


# ============================================================
# 2. 设置导出文件夹
# ============================================================

output_dir <- input_folder

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 3. 设置 marker 名称
# ============================================================
# 要算 RGT 就写 "RGT"
# 要算 LGT 就写 "LGT"

marker_name <- "RGT"


# ============================================================
# 4. 参数设置
# ============================================================

force_threshold <- 20
forceplate_index <- 1

# 如果 Mokka 第一帧不是 1，比如是 307，就写 307
# 如果想自动识别，就保持 NA
first_mokka_frame_manual <- NA

# 是否递归处理子文件夹
recursive_search <- TRUE

# 是否自动排除 static 文件
exclude_static_files <- TRUE


# ============================================================
# 5. 自动识别起始帧函数
# ============================================================

get_first_mokka_frame_auto <- function(c3d_object) {
  
  get_by_path <- function(x, path) {
    current <- x
    for (p in path) {
      if (is.list(current) && p %in% names(current)) {
        current <- current[[p]]
      } else {
        return(NULL)
      }
    }
    return(current)
  }
  
  extract_num <- function(x) {
    if (is.null(x)) return(NA_real_)
    
    if (is.list(x)) {
      if ("data" %in% names(x)) {
        x <- x$data
      } else if ("value" %in% names(x)) {
        x <- x$value
      } else if ("values" %in% names(x)) {
        x <- x$values
      } else {
        x <- unlist(x)
      }
    }
    
    x_num <- suppressWarnings(as.numeric(unlist(x)))
    x_num <- x_num[is.finite(x_num)]
    
    if (length(x_num) == 0) return(NA_real_)
    return(x_num[1])
  }
  
  possible_paths <- list(
    c("header", "first_frame"),
    c("header", "firstFrame"),
    c("header", "frame_first"),
    c("header", "point", "first_frame"),
    c("header", "points", "first_frame"),
    c("parameters", "POINT", "FIRST_FRAME"),
    c("parameters", "POINT", "FIRST_FRAME", "data"),
    c("parameters", "POINT", "FIRST_FRAME", "value"),
    c("parameters", "POINT", "FIRST_FRAME", "values")
  )
  
  for (path in possible_paths) {
    value <- get_by_path(c3d_object, path)
    value_num <- extract_num(value)
    
    if (!is.na(value_num) && value_num > 0 && value_num < 1000000) {
      return(as.integer(value_num))
    }
  }
  
  return(NA_integer_)
}


# ============================================================
# 6. 定义结果表列名
# ============================================================

result_cols <- c(
  "Status",
  "Error_Message",
  "C3D_File",
  "Full_C3D_Path",
  "Parent_Folder",
  
  "First_Mokka_Frame",
  
  "Forceplate_Index",
  "Force_Threshold_N",
  "Analog_Per_Frame",
  
  "Pre_Onset_Force_Index",
  "Pre_Onset_R_Row",
  "Pre_Onset_Mokka_Frame",
  "Pre_Onset_Force_Raw",
  "Pre_Onset_Force_Abs",
  
  "Onset_Force_Index",
  "Onset_R_Row",
  "Onset_Mokka_Frame",
  "Onset_Force_Raw",
  "Onset_Force_Abs",
  
  "Offset_Force_Index",
  "Offset_R_Row",
  "Offset_Mokka_Frame",
  "Offset_Force_Raw",
  "Offset_Force_Abs",
  
  "Peak_Force_Index",
  "Peak_R_Row",
  "Peak_Mokka_Frame",
  "Peak_Force_Raw",
  "Peak_Force_Abs",
  
  "Marker",
  
  "Marker_X_Onset",
  "Marker_Y_Onset",
  "Marker_Z_Onset",
  
  "Marker_X_PeakForce",
  "Marker_Y_PeakForce",
  "Marker_Z_PeakForce",
  
  "PeakForce_Z_Displacement",
  "PeakForce_Z_Displacement_Abs",
  
  "Marker_Frame_Difference_Onset_to_PeakForce",
  "Mokka_Frame_Difference_Onset_to_PeakForce",
  
  "PeakForceAbs_Div_PeakForceZDisplacementAbs_N_per_mm",
  "PeakForceAbs_Div_PeakForceZDisplacementAbs_N_per_m",
  "PeakForceAbs_Div_PeakForceZDisplacementAbs_kN_per_m"
)


make_empty_result <- function(c3d_file, status = "ERROR", error_message = NA_character_) {
  
  x <- as.list(rep(NA, length(result_cols)))
  names(x) <- result_cols
  
  x$Status <- status
  x$Error_Message <- error_message
  x$C3D_File <- basename(c3d_file)
  x$Full_C3D_Path <- c3d_file
  x$Parent_Folder <- dirname(c3d_file)
  x$Marker <- marker_name
  x$Forceplate_Index <- forceplate_index
  x$Force_Threshold_N <- force_threshold
  
  return(as.data.frame(x, check.names = FALSE))
}


# ============================================================
# 7. 单个 C3D 文件处理函数
# ============================================================

process_one_c3d <- function(c3d_file) {
  
  cat("\n============================================================\n")
  cat("正在处理：\n")
  cat(c3d_file, "\n")
  cat("============================================================\n")
  
  tryCatch({
    
    # ---------- 读取 C3D ----------
    d <- c3d_read(c3d_file)
    
    
    # ---------- 起始帧 ----------
    if (!is.na(first_mokka_frame_manual)) {
      first_mokka_frame <- first_mokka_frame_manual
    } else {
      first_mokka_frame <- get_first_mokka_frame_auto(d)
      
      if (is.na(first_mokka_frame)) {
        first_mokka_frame <- 1
      }
    }
    
    
    # ---------- marker 数据 ----------
    marker_data <- c3d_data(d, format = "wide")
    
    n_marker_frames <- nrow(marker_data)
    
    if (is.null(marker_data) || n_marker_frames == 0) {
      stop("marker_data 为空，无法分析。")
    }
    
    
    # ---------- forceplatform 数据 ----------
    if (is.null(d$forceplatform) || length(d$forceplatform) == 0) {
      stop("该 C3D 没有读取到 forceplatform 数据。")
    }
    
    if (forceplate_index > length(d$forceplatform)) {
      stop("forceplate_index 超过该 C3D 中的力平台数量。")
    }
    
    fp <- d$forceplatform[[forceplate_index]]
    
    if (is.null(fp$forces)) {
      stop("当前 force platform 中没有 forces 数据。")
    }
    
    force_matrix <- fp$forces
    
    if (ncol(force_matrix) < 3) {
      stop("force_matrix 少于 3 列，无法提取第 3 列 Fz。")
    }
    
    # 第 3 列一般为垂直方向力 Fz
    force_raw_highfreq <- as.numeric(force_matrix[, 3])
    force_abs_highfreq <- abs(force_raw_highfreq)
    
    n_force_frames <- length(force_abs_highfreq)
    
    if (all(is.na(force_abs_highfreq))) {
      stop("力数据全部为 NA，无法分析。")
    }
    
    
    # ---------- analog_per_frame ----------
    analog_per_frame <- NA
    
    if (!is.null(d$header$analogperframe)) {
      analog_per_frame <- suppressWarnings(as.numeric(d$header$analogperframe))
    }
    
    if (is.na(analog_per_frame) || analog_per_frame <= 0) {
      analog_per_frame <- round(n_force_frames / n_marker_frames)
    }
    
    analog_per_frame <- as.integer(round(analog_per_frame))
    
    if (is.na(analog_per_frame) || analog_per_frame <= 0) {
      stop("analog_per_frame 无法识别。")
    }
    
    
    # ========================================================
    # 找最大垂直地面反作用力
    # ========================================================
    
    peak_force_index <- which.max(force_abs_highfreq)
    
    peak_force_raw <- force_raw_highfreq[peak_force_index]
    peak_force_abs <- force_abs_highfreq[peak_force_index]
    
    peak_row_index <- ceiling(peak_force_index / analog_per_frame)
    peak_row_index <- max(1, min(peak_row_index, n_marker_frames))
    
    peak_mokka_frame <- first_mokka_frame + peak_row_index - 1
    
    
    # ---------- 找力开始出现帧 ----------
    force_present <- force_abs_highfreq > force_threshold
    
    if (!force_present[peak_force_index]) {
      stop("力峰值帧没有超过阈值，请检查 force_threshold。")
    }
    
    onset_force_index <- peak_force_index
    
    while (onset_force_index > 1 && force_present[onset_force_index - 1]) {
      onset_force_index <- onset_force_index - 1
    }
    
    onset_force_raw <- force_raw_highfreq[onset_force_index]
    onset_force_abs <- force_abs_highfreq[onset_force_index]
    
    onset_row_index <- ceiling(onset_force_index / analog_per_frame)
    onset_row_index <- max(1, min(onset_row_index, n_marker_frames))
    
    onset_mokka_frame <- first_mokka_frame + onset_row_index - 1
    
    
    # ---------- 找力结束帧 ----------
    offset_force_index <- peak_force_index
    
    while (offset_force_index < n_force_frames && force_present[offset_force_index + 1]) {
      offset_force_index <- offset_force_index + 1
    }
    
    offset_force_raw <- force_raw_highfreq[offset_force_index]
    offset_force_abs <- force_abs_highfreq[offset_force_index]
    
    offset_row_index <- ceiling(offset_force_index / analog_per_frame)
    offset_row_index <- max(1, min(offset_row_index, n_marker_frames))
    
    offset_mokka_frame <- first_mokka_frame + offset_row_index - 1
    
    
    # ---------- 找力开始出现前一帧 ----------
    pre_onset_force_index <- onset_force_index - 1
    
    if (pre_onset_force_index >= 1) {
      
      pre_onset_force_raw <- force_raw_highfreq[pre_onset_force_index]
      pre_onset_force_abs <- force_abs_highfreq[pre_onset_force_index]
      
      pre_onset_row_index <- ceiling(pre_onset_force_index / analog_per_frame)
      pre_onset_row_index <- max(1, min(pre_onset_row_index, n_marker_frames))
      
      pre_onset_mokka_frame <- first_mokka_frame + pre_onset_row_index - 1
      
    } else {
      
      pre_onset_force_raw <- NA_real_
      pre_onset_force_abs <- NA_real_
      pre_onset_row_index <- NA_real_
      pre_onset_mokka_frame <- NA_real_
    }
    
    
    # ---------- 找 marker 的 X/Y/Z 坐标列 ----------
    x_col <- paste0(marker_name, "_x")
    y_col <- paste0(marker_name, "_y")
    z_col <- paste0(marker_name, "_z")
    
    needed_cols <- c(x_col, y_col, z_col)
    
    if (!all(needed_cols %in% names(marker_data))) {
      
      possible_cols <- names(marker_data)[
        grepl(
          paste0("^", marker_name, "_"),
          names(marker_data),
          ignore.case = TRUE
        )
      ]
      
      if (length(possible_cols) < 3) {
        stop(paste0("没有找到 marker ", marker_name, " 的完整 X/Y/Z 坐标。"))
      }
      
      x_col <- possible_cols[grepl("_x$", possible_cols, ignore.case = TRUE)][1]
      y_col <- possible_cols[grepl("_y$", possible_cols, ignore.case = TRUE)][1]
      z_col <- possible_cols[grepl("_z$", possible_cols, ignore.case = TRUE)][1]
      
      needed_cols <- c(x_col, y_col, z_col)
    }
    
    
    # ========================================================
    # 提取最大力时刻对应的 marker 位移
    # ========================================================
    
    coord_onset <- marker_data[onset_row_index, needed_cols, drop = FALSE]
    
    x_onset <- as.numeric(coord_onset[[x_col]])
    y_onset <- as.numeric(coord_onset[[y_col]])
    z_onset <- as.numeric(coord_onset[[z_col]])
    
    coord_peak <- marker_data[peak_row_index, needed_cols, drop = FALSE]
    
    x_peak <- as.numeric(coord_peak[[x_col]])
    y_peak <- as.numeric(coord_peak[[y_col]])
    z_peak <- as.numeric(coord_peak[[z_col]])
    
    peakforce_z_displacement <- z_peak - z_onset
    peakforce_z_displacement_abs <- abs(peakforce_z_displacement)
    
    marker_frame_difference_onset_to_peakforce <- peak_row_index - onset_row_index
    mokka_frame_difference_onset_to_peakforce <- peak_mokka_frame - onset_mokka_frame
    
    
    # ========================================================
    # 计算刚度
    # ========================================================
    # 刚度 = 最大垂直地面反作用力绝对值 / 最大力时刻对应的Z轴位移绝对值
    
    if (is.na(peakforce_z_displacement_abs) || peakforce_z_displacement_abs == 0) {
      
      peak_force_abs_div_peak_zdisp_abs_N_per_mm <- NA_real_
      peak_force_abs_div_peak_zdisp_abs_N_per_m <- NA_real_
      peak_force_abs_div_peak_zdisp_abs_kN_per_m <- NA_real_
      
    } else {
      
      peak_force_abs_div_peak_zdisp_abs_N_per_mm <- 
        peak_force_abs / peakforce_z_displacement_abs
      
      peak_force_abs_div_peak_zdisp_abs_N_per_m <- 
        peak_force_abs_div_peak_zdisp_abs_N_per_mm * 1000
      
      peak_force_abs_div_peak_zdisp_abs_kN_per_m <- 
        peak_force_abs_div_peak_zdisp_abs_N_per_mm
    }
    
    
    # ---------- 生成该文件结果 ----------
    result <- data.frame(
      Status = "OK",
      Error_Message = NA_character_,
      C3D_File = basename(c3d_file),
      Full_C3D_Path = c3d_file,
      Parent_Folder = dirname(c3d_file),
      
      First_Mokka_Frame = first_mokka_frame,
      
      Forceplate_Index = forceplate_index,
      Force_Threshold_N = force_threshold,
      Analog_Per_Frame = analog_per_frame,
      
      Pre_Onset_Force_Index = pre_onset_force_index,
      Pre_Onset_R_Row = pre_onset_row_index,
      Pre_Onset_Mokka_Frame = pre_onset_mokka_frame,
      Pre_Onset_Force_Raw = pre_onset_force_raw,
      Pre_Onset_Force_Abs = pre_onset_force_abs,
      
      Onset_Force_Index = onset_force_index,
      Onset_R_Row = onset_row_index,
      Onset_Mokka_Frame = onset_mokka_frame,
      Onset_Force_Raw = onset_force_raw,
      Onset_Force_Abs = onset_force_abs,
      
      Offset_Force_Index = offset_force_index,
      Offset_R_Row = offset_row_index,
      Offset_Mokka_Frame = offset_mokka_frame,
      Offset_Force_Raw = offset_force_raw,
      Offset_Force_Abs = offset_force_abs,
      
      Peak_Force_Index = peak_force_index,
      Peak_R_Row = peak_row_index,
      Peak_Mokka_Frame = peak_mokka_frame,
      Peak_Force_Raw = peak_force_raw,
      Peak_Force_Abs = peak_force_abs,
      
      Marker = marker_name,
      
      Marker_X_Onset = x_onset,
      Marker_Y_Onset = y_onset,
      Marker_Z_Onset = z_onset,
      
      Marker_X_PeakForce = x_peak,
      Marker_Y_PeakForce = y_peak,
      Marker_Z_PeakForce = z_peak,
      
      PeakForce_Z_Displacement = peakforce_z_displacement,
      PeakForce_Z_Displacement_Abs = peakforce_z_displacement_abs,
      
      Marker_Frame_Difference_Onset_to_PeakForce = 
        marker_frame_difference_onset_to_peakforce,
      
      Mokka_Frame_Difference_Onset_to_PeakForce = 
        mokka_frame_difference_onset_to_peakforce,
      
      PeakForceAbs_Div_PeakForceZDisplacementAbs_N_per_mm = 
        peak_force_abs_div_peak_zdisp_abs_N_per_mm,
      
      PeakForceAbs_Div_PeakForceZDisplacementAbs_N_per_m = 
        peak_force_abs_div_peak_zdisp_abs_N_per_m,
      
      PeakForceAbs_Div_PeakForceZDisplacementAbs_kN_per_m = 
        peak_force_abs_div_peak_zdisp_abs_kN_per_m,
      
      check.names = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    
    cat("处理失败：", basename(c3d_file), "\n")
    cat("错误原因：", conditionMessage(e), "\n")
    
    return(make_empty_result(
      c3d_file = c3d_file,
      status = "ERROR",
      error_message = conditionMessage(e)
    ))
  })
}


# ============================================================
# 8. 批量寻找 C3D 文件
# ============================================================

c3d_files <- list.files(
  path = input_folder,
  pattern = "\\.c3d$",
  full.names = TRUE,
  recursive = recursive_search,
  ignore.case = TRUE
)

if (exclude_static_files) {
  c3d_files <- c3d_files[
    !grepl("static", basename(c3d_files), ignore.case = TRUE)
  ]
}

if (length(c3d_files) == 0) {
  stop("该文件夹下没有找到需要处理的 C3D 文件。")
}

cat("\n============================================================\n")
cat("共找到 C3D 文件数量：", length(c3d_files), "\n")
cat("Marker：", marker_name, "\n")
cat("输入文件夹：", input_folder, "\n")
cat("是否搜索子文件夹：", recursive_search, "\n")
cat("是否排除 static 文件：", exclude_static_files, "\n")
cat("============================================================\n\n")

print(c3d_files)


# ============================================================
# 9. 批量处理所有 C3D
# ============================================================

all_results_list <- lapply(c3d_files, process_one_c3d)

all_results <- do.call(rbind, all_results_list)


# ============================================================
# 10. 导出一个总表
# ============================================================

out_summary <- file.path(
  output_dir,
  paste0(
    "POGO_CMJ_", marker_name,
    "_batch_peakForce_to_peakForceDisplacement_results.csv"
  )
)

write.csv(all_results, out_summary, row.names = FALSE)

cat("\n============================================================\n")
cat("POGO CMJ 分析完成！\n\n")
cat("成功数量：", sum(all_results$Status == "OK"), "\n")
cat("失败数量：", sum(all_results$Status == "ERROR"), "\n\n")
cat("总表已导出到：\n")
cat(out_summary, "\n")
cat("============================================================\n")

shell.exec(output_dir)