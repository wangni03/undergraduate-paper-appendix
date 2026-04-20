#' 深入调试模型结构
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param max_depth 最大递归深度
deep_debug_model <- function(model_obj, model_name, max_depth = 3) {
  cat(sprintf("\n=== 深度调试模型: %s ===\n", model_name))
  cat(sprintf("对象类型: %s\n", paste(class(model_obj), collapse=", ")))
  
  if (is.list(model_obj)) {
    cat(sprintf("列表长度: %d\n", length(model_obj)))
    cat("顶层字段名: ", paste(names(model_obj), collapse=", "), "\n")
    
    # 检查是否有fit字段
    if ("fit" %in% names(model_obj)) {
      cat("\n--- fit字段结构 ---\n")
      fit_obj <- model_obj$fit
      cat(sprintf("fit对象类型: %s\n", paste(class(fit_obj), collapse=", ")))
      
      if (is.list(fit_obj)) {
        cat("fit字段名: ", paste(names(fit_obj), collapse=", "), "\n")
      }
    }
    
    # 检查是否有model字段
    if ("model" %in% names(model_obj)) {
      cat("\n--- model字段结构 ---\n")
      model_field <- model_obj$model
      cat(sprintf("model对象类型: %s\n", paste(class(model_field), collapse=", ")))
    }
    
    # 检查是否有features或featurenames
    feature_fields <- c("features", "featurenames", "featureNames", "geneNames", "var.names")
    for (field in feature_fields) {
      if (field %in% names(model_obj)) {
        cat(sprintf("\n--- %s字段 ---\n", field))
        field_value <- model_obj[[field]]
        if (is.character(field_value)) {
          cat(sprintf("类型: character, 长度: %d\n", length(field_value)))
          cat("示例: ", head(field_value, 5), "\n")
        } else {
          cat(sprintf("类型: %s\n", class(field_value)))
        }
      }
    }
    
    # 检查是否有importance字段
    if ("importance" %in% names(model_obj)) {
      cat("\n--- importance字段 ---\n")
      imp <- model_obj$importance
      cat(sprintf("类型: %s\n", class(imp)))
      if (!is.null(names(imp))) {
        cat("有名称: 是\n")
        cat("长度: ", length(imp), "\n")
        cat("示例: ", head(names(imp), 5), "\n")
      }
    }
    
    # 检查是否有coefficients字段
    if ("coefficients" %in% names(model_obj)) {
      cat("\n--- coefficients字段 ---\n")
      coefs <- model_obj$coefficients
      cat(sprintf("类型: %s\n", class(coefs)))
      if (!is.null(names(coefs))) {
        cat("有名称: 是\n")
        cat("非零系数数量: ", sum(coefs != 0, na.rm = TRUE), "\n")
        cat("示例: ", head(names(coefs)[which(coefs != 0)], 5), "\n")
      }
    }
    
    # 检查是否有x或data字段（训练数据）
    data_fields <- c("x", "X", "data", "dataX", "train", "trainX")
    for (field in data_fields) {
      if (field %in% names(model_obj)) {
        cat(sprintf("\n--- %s字段（训练数据）---\n", field))
        data_obj <- model_obj[[field]]
        cat(sprintf("类型: %s\n", class(data_obj)))
        
        if (is.matrix(data_obj)) {
          cat(sprintf("维度: %s\n", paste(dim(data_obj), collapse=" x ")))
          if (!is.null(colnames(data_obj))) {
            cat("有列名: 是\n")
            cat("列数: ", ncol(data_obj), "\n")
            cat("列名示例: ", head(colnames(data_obj), 5), "\n")
          }
        } else if (is.data.frame(data_obj)) {
          cat(sprintf("维度: %s\n", paste(dim(data_obj), collapse=" x ")))
          cat("列数: ", length(data_obj), "\n")
          cat("列名示例: ", head(names(data_obj), 5), "\n")
        }
      }
    }
    
    # 递归检查子对象（限制深度）
    if (max_depth > 0) {
      cat(sprintf("\n--- 递归检查子对象（深度限制: %d）---\n", max_depth))
      for (i in seq_along(model_obj)) {
        elem_name <- names(model_obj)[i]
        if (is.null(elem_name) || elem_name == "") {
          elem_name <- paste0("[[", i, "]]")
        }
        
        elem <- model_obj[[i]]
        
        # 只检查列表对象
        if (is.list(elem) && length(elem) > 0 && max_depth > 1) {
          # 检查是否可能是模型对象
          if (any(sapply(elem, function(x) inherits(x, c("rfsrc", "coxph", "cv.glmnet", "CoxBoost", "plsRcox", "survivalsvm"))))) {
            cat(sprintf("\n%s 包含模型对象:\n", elem_name))
            for (j in seq_along(elem)) {
              sub_elem <- elem[[j]]
              if (inherits(sub_elem, c("rfsrc", "coxph", "cv.glmnet", "CoxBoost", "plsRcox", "survivalsvm"))) {
                cat(sprintf("  %s: %s\n", names(elem)[j], class(sub_elem)[1]))
              }
            }
          }
        }
      }
    }
  }
  
  cat("\n=== 调试结束 ===\n")
}

#' 从ML.Dev.Prog.Sig结果中提取每个模型的基因名称
#'
#' @param res ML.Dev.Prog.Sig函数返回的结果对象
#' @param verbose 是否显示调试信息，默认为FALSE
#' @param debug_missing 是否对未能提取基因的模型进行深度调试，默认为TRUE
#' @return 包含每个模型及其对应基因名称的列表
#' @export
#'
#' @examples
get_genes_by_model <- function(res, verbose = FALSE, debug_missing = TRUE) {
  if (is.null(res$ml.res)) {
    stop("结果对象中没有找到ml.res")
  }
  
  # 调试：打印总模型数
  total_models <- length(res$ml.res)
  if (verbose) {
    cat(sprintf("结果中共有 %d 个模型\n", total_models))
  }
  
  # 初始化结果列表
  genes_by_model <- list()
  missing_models <- c()
  
  # 遍历所有模型
  for (model_name in names(res$ml.res)) {
    model_obj <- res$ml.res[[model_name]]
    
    # 对特定模型进行深入调试
    if (model_name %in% c("StepCox[both] + plsRcox", "StepCox[forward] + plsRcox", 
                         "StepCox[backward] + plsRcox", "Lasso + plsRcox",
                         "CoxBoost + plsRcox", "RSF + plsRcox")) {
      cat(sprintf("\n=== 深入分析模型: %s ===\n", model_name))
      genes <- extract_genes_from_model_with_debug(model_obj, model_name, verbose = TRUE)
    } else {
      genes <- extract_genes_from_model(model_obj, model_name, verbose)
    }
    
    if (is.null(genes) || length(genes) == 0) {
      missing_models <- c(missing_models, model_name)
      genes_by_model[[model_name]] <- character(0)  # 空向量而不是NULL
      
      # 如果需要调试缺失的模型
      if (debug_missing && model_name %in% missing_models) {
        deep_debug_model(model_obj, model_name, max_depth = 3)
      }
    } else {
      genes_by_model[[model_name]] <- genes
    }
  }
  
  # 打印警告信息
  if (length(missing_models) > 0 && verbose) {
    cat(sprintf("\n警告：%d 个模型未能提取基因: %s\n", 
                length(missing_models),
                paste(head(missing_models, 5), collapse=", ")))
    if (length(missing_models) > 5) {
      cat(sprintf("... 以及另外 %d 个模型\n", length(missing_models) - 5))
    }
  }
  
  if (verbose) {
    cat(sprintf("\n成功处理了 %d/%d 个模型\n", 
                length(genes_by_model) - length(missing_models), 
                total_models))
  }
  
  # 返回结果
  return(genes_by_model)
}

#' 从单个模型对象中提取基因名称（带调试）
#'
#' @param model_obj 单个模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息，默认为FALSE
#' @return 该模型选择的基因名称向量
extract_genes_from_model_with_debug <- function(model_obj, model_name, verbose = FALSE) {
  
  # 调试：打印处理的模型名称
  cat(sprintf("\n=== 深度分析模型: %s ===\n", model_name))
  cat(sprintf("对象类型: %s\n", paste(class(model_obj), collapse=", ")))
  
  if (is.list(model_obj)) {
    cat(sprintf("列表长度: %d\n", length(model_obj)))
    cat("顶层字段名: ", paste(names(model_obj), collapse=", "), "\n")
    
    # 遍历所有元素
    for (i in seq_along(model_obj)) {
      elem_name <- names(model_obj)[i]
      if (is.null(elem_name) || elem_name == "") {
        elem_name <- paste0("[[", i, "]]")
      }
      
      elem <- model_obj[[i]]
      cat(sprintf("\n元素 %s:\n", elem_name))
      cat(sprintf("  类型: %s\n", paste(class(elem), collapse=", ")))
      
      if (is.list(elem)) {
        cat(sprintf("  列表长度: %d\n", length(elem)))
        if (length(elem) > 0) {
          cat(sprintf("  字段名: %s\n", paste(names(elem), collapse=", ")))
        }
        
        # 检查是否是模型对象
        if (inherits(elem, "coxph")) {
          cat("  -> 发现coxph对象\n")
          if (!is.null(elem$coefficients)) {
            cat(sprintf("     系数数量: %d\n", length(elem$coefficients)))
            non_zero <- which(elem$coefficients != 0)
            if (length(non_zero) > 0) {
              cat("     非零系数基因: ", names(elem$coefficients)[non_zero], "\n")
            } else {
              cat("     所有系数都为零\n")
            }
          }
        } else if (inherits(elem, "plsRcox") || inherits(elem, "plsRcoxmodel")) {
          cat("  -> 发现plsRcox对象\n")
          # 检查plsRcox的字段
          cat("     字段: ", paste(names(elem), collapse=", "), "\n")
          
          # 特别检查dataX字段
          if ("dataX" %in% names(elem)) {
            dataX <- elem$dataX
            if (is.data.frame(dataX)) {
              cat(sprintf("     dataX列数: %d\n", length(dataX)))
              cat("     列名: ", names(dataX), "\n")
            } else if (is.matrix(dataX) && !is.null(colnames(dataX))) {
              cat(sprintf("     dataX列数: %d\n", ncol(dataX)))
              cat("     列名: ", colnames(dataX), "\n")
            }
          }
          
          if ("ExpliX" %in% names(elem)) {
            expliX <- elem$ExpliX
            if (is.matrix(expliX) && !is.null(colnames(expliX))) {
              cat(sprintf("     ExpliX列数: %d\n", ncol(expliX)))
              cat("     列名: ", head(colnames(expliX), 5), "\n")
            }
          }
          
          if ("Xnames" %in% names(elem)) {
            cat(sprintf("     Xnames长度: %d\n", length(elem$Xnames)))
            cat("     示例: ", head(elem$Xnames, 5), "\n")
          }
          if ("var.names" %in% names(elem)) {
            cat(sprintf("     var.names长度: %d\n", length(elem$var.names)))
            cat("     示例: ", head(elem$var.names, 5), "\n")
          }
        } else if (inherits(elem, "cv.glmnet")) {
          cat("  -> 发现cv.glmnet对象\n")
          tryCatch({
            coef_mat <- as.matrix(coef(elem, s = "lambda.min"))
            non_zero <- which(coef_mat != 0)
            cat(sprintf("     总系数: %d, 非零系数: %d\n", length(coef_mat), length(non_zero)))
            if (length(non_zero) > 0) {
              genes <- rownames(coef_mat)[non_zero]
              genes <- genes[!grepl("Intercept", genes)]
              cat("     非零系数基因: ", head(genes, 5), "\n")
            }
          }, error = function(e) cat("     无法提取系数\n"))
        } else if (inherits(elem, "CoxBoost")) {
          cat("  -> 发现CoxBoost对象\n")
          tryCatch({
            coef_mat <- coef(elem)
            non_zero <- which(coef_mat != 0)
            cat(sprintf("     总系数: %d, 非零系数: %d\n", length(coef_mat), length(non_zero)))
            if (length(non_zero) > 0) {
              cat("     非零系数基因: ", head(names(coef_mat)[non_zero], 5), "\n")
            }
          }, error = function(e) cat("     无法提取系数\n"))
        } else if (inherits(elem, "rfsrc")) {
          cat("  -> 发现rfsrc对象\n")
          if (!is.null(elem$importance)) {
            cat(sprintf("     重要性变量数量: %d\n", length(elem$importance)))
            cat("     示例变量: ", head(names(elem$importance), 5), "\n")
          }
        }
      }
    }
  }
  
  # 现在尝试提取基因
  cat("\n=== 尝试提取基因 ===\n")
  
  # 对于plsRcoxmodel对象，直接从dataX提取
  if (inherits(model_obj, "plsRcoxmodel")) {
    cat("对象是plsRcoxmodel类型，尝试从dataX字段提取基因\n")
    if (!is.null(model_obj$dataX)) {
      if (is.data.frame(model_obj$dataX)) {
        genes <- names(model_obj$dataX)
        if (length(genes) > 0) {
          cat(sprintf("从dataX字段提取了 %d 个基因: %s\n", 
                      length(genes), 
                      paste(genes, collapse=", ")))
          return(genes)
        }
      } else if (is.matrix(model_obj$dataX) && !is.null(colnames(model_obj$dataX))) {
        genes <- colnames(model_obj$dataX)
        if (length(genes) > 0) {
          cat(sprintf("从dataX矩阵提取了 %d 个基因: %s\n", 
                      length(genes), 
                      paste(genes, collapse=", ")))
          return(genes)
        }
      }
    }
  }
  
  # 首先检查是否是StepCox + plsRcox
  if (grepl("StepCox\\[", model_name) && grepl(" \\+ plsRcox", model_name)) {
    cat("检测到StepCox + plsRcox组合\n")
    
    # 方法1: 尝试从plsRcox部分提取原始基因
    plsRcox_genes <- find_and_extract_plsRcox_genes(model_obj, TRUE)
    if (!is.null(plsRcox_genes) && length(plsRcox_genes) > 0) {
      cat(sprintf("从plsRcox部分提取了 %d 个基因\n", length(plsRcox_genes)))
      return(plsRcox_genes)
    }
    
    # 方法2: 尝试直接查找coxph对象
    coxph_genes <- find_and_extract_coxph_genes(model_obj, TRUE)
    if (!is.null(coxph_genes) && length(coxph_genes) > 0) {
      cat(sprintf("从coxph对象提取了 %d 个基因: %s\n", 
                  length(coxph_genes), 
                  paste(head(coxph_genes, 5), collapse=", ")))
      return(coxph_genes)
    }
    
    # 方法3: 尝试从训练数据中提取
    train_data_genes <- find_and_extract_training_genes(model_obj, TRUE)
    if (!is.null(train_data_genes) && length(train_data_genes) > 0) {
      cat(sprintf("从训练数据提取了 %d 个基因\n", length(train_data_genes)))
      return(train_data_genes)
    }
  }
  
  # 对于其他包含plsRcox的组合
  if (grepl(" \\+ plsRcox", model_name)) {
    cat("检测到其他plsRcox组合\n")
    
    # 提取主组件名称
    main_component <- gsub(" \\+ plsRcox", "", model_name)
    cat(sprintf("主组件: %s\n", main_component))
    
    # 根据主组件类型尝试提取
    if (grepl("Lasso", main_component)) {
      main_genes <- try_extract_lasso_component(model_obj, TRUE)
    } else if (grepl("CoxBoost", main_component)) {
      main_genes <- try_extract_coxboost_component(model_obj, TRUE)
    } else if (grepl("RSF", main_component)) {
      main_genes <- try_extract_rsf_component(model_obj, TRUE)
    } else {
      main_genes <- NULL
    }
    
    if (!is.null(main_genes) && length(main_genes) > 0) {
      cat(sprintf("从主组件提取了 %d 个基因\n", length(main_genes)))
      return(main_genes)
    }
    
    # 如果主组件没有提取到，尝试从plsRcox部分提取
    plsRcox_genes <- find_and_extract_plsRcox_genes(model_obj, TRUE)
    if (!is.null(plsRcox_genes) && length(plsRcox_genes) > 0) {
      cat(sprintf("从plsRcox部分提取了 %d 个特征\n", length(plsRcox_genes)))
      return(plsRcox_genes)
    }
  }
  
  cat("=== 无法提取基因 ===\n")
  return(NULL)
}


#' 查找并提取coxph对象的基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
find_and_extract_coxph_genes <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("查找coxph对象\n")
  
  # 查找coxph对象
  coxph_obj <- find_object_by_class(model_obj, "coxph")
  if (!is.null(coxph_obj)) {
    if (verbose) cat("找到coxph对象\n")
    genes <- names(which(!is.na(coef(coxph_obj)) & coef(coxph_obj) != 0))
    if (length(genes) > 0) {
      return(genes)
    }
  }
  
  # 检查是否有fit字段且是coxph
  if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "coxph")) {
    if (verbose) cat("在fit字段找到coxph对象\n")
    genes <- names(which(!is.na(coef(model_obj$fit)) & coef(model_obj$fit) != 0))
    if (length(genes) > 0) {
      return(genes)
    }
  }
  
  # 检查是否有model字段且是coxph
  if (is.list(model_obj) && "model" %in% names(model_obj) && inherits(model_obj$model, "coxph")) {
    if (verbose) cat("在model字段找到coxph对象\n")
    genes <- names(which(!is.na(coef(model_obj$model)) & coef(model_obj$model) != 0))
    if (length(genes) > 0) {
      return(genes)
    }
  }
  
  return(NULL)
}

#' 查找并提取训练数据中的基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
find_and_extract_training_genes <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("查找训练数据中的基因\n")
  
  if (is.list(model_obj)) {
    # 检查训练数据字段
    train_fields <- c("x", "X", "data", "train.x", "trainX", "xdata")
    for (field in train_fields) {
      if (field %in% names(model_obj)) {
        data_obj <- model_obj[[field]]
        if (is.matrix(data_obj) && !is.null(colnames(data_obj))) {
          if (verbose) cat(sprintf("从 %s 字段提取基因\n", field))
          return(colnames(data_obj))
        } else if (is.data.frame(data_obj) && length(data_obj) > 0) {
          if (verbose) cat(sprintf("从 %s 字段提取基因\n", field))
          return(names(data_obj))
        }
      }
    }
  }
  
  return(NULL)
}

#' 查找并提取plsRcox对象的基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
find_and_extract_plsRcox_genes <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("查找plsRcox对象\n")
  
  # 查找plsRcox对象
  plsRcox_obj <- find_object_by_class(model_obj, "plsRcox")
  if (!is.null(plsRcox_obj)) {
    if (verbose) cat("找到plsRcox对象\n")
    genes <- extract_plsRcox_genes_safely(plsRcox_obj, verbose)
    return(genes)
  }
  
  # 查找plsRcoxmodel对象
  plsRcoxmodel_obj <- find_object_by_class(model_obj, "plsRcoxmodel")
  if (!is.null(plsRcoxmodel_obj)) {
    if (verbose) cat("找到plsRcoxmodel对象\n")
    genes <- extract_genes_from_plsRcoxmodel(plsRcoxmodel_obj, "plsRcoxmodel", verbose)
    return(genes)
  }
  
  # 检查是否有fit字段且是plsRcox
  if (is.list(model_obj) && "fit" %in% names(model_obj)) {
    if (inherits(model_obj$fit, "plsRcox")) {
      if (verbose) cat("在fit字段找到plsRcox对象\n")
      genes <- extract_plsRcox_genes_safely(model_obj$fit, verbose)
      return(genes)
    } else if (inherits(model_obj$fit, "plsRcoxmodel")) {
      if (verbose) cat("在fit字段找到plsRcoxmodel对象\n")
      genes <- extract_genes_from_plsRcoxmodel(model_obj$fit, "plsRcoxmodel", verbose)
      return(genes)
    }
  }
  
  return(NULL)
}


#' 从单个模型对象中提取基因名称
#'
#' @param model_obj 单个模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息，默认为FALSE
#' @return 该模型选择的基因名称向量
extract_genes_from_model <- function(model_obj, model_name, verbose = FALSE) {
  
  # 对于包含plsRcox的组合模型，使用专门的函数
  if (model_name %in% c("StepCox[both] + plsRcox", "StepCox[forward] + plsRcox", 
                       "StepCox[backward] + plsRcox", "Lasso + plsRcox",
                       "CoxBoost + plsRcox", "RSF + plsRcox")) {
    if (verbose) cat(sprintf("使用专门函数处理 %s\n", model_name))
    return(extract_genes_from_model_with_debug(model_obj, model_name, verbose))
  }
  
  # 调试：打印处理的模型名称
  if (verbose) {
    cat(sprintf("\n处理模型: %s\n", model_name))
    cat(sprintf("对象类型: %s\n", paste(class(model_obj), collapse=", ")))
  }
  
  # 首先检查是否是组合模型
  if (grepl(" \\+ ", model_name) || grepl("\\[both\\] \\+ ", model_name)) {
    if (verbose) {
      cat("识别为组合模型\n")
    }
    
    # 对于需要特殊处理的模型
    special_models <- c("RSF + SuperPC", "RSF + survival-SVM", "StepCox[both] + SuperPC",
                       "StepCox[backward] + SuperPC", "StepCox[forward] + SuperPC",
                       "CoxBoost + SuperPC", "Lasso + SuperPC",
                       "StepCox[both] + RSF", "StepCox[backward] + RSF", 
                       "StepCox[forward] + RSF", "Lasso + RSF", 
                       "Lasso + plsRcox", "CoxBoost + plsRcox", "RSF + plsRcox",
                       "StepCox[both] + plsRcox", "StepCox[forward] + plsRcox",
                       "StepCox[backward] + plsRcox")
    
    if (model_name %in% special_models) {
      if (verbose) cat("这是需要特殊处理的模型\n")
      genes <- extract_genes_from_special_ensemble(model_obj, model_name, verbose)
      if (!is.null(genes) && length(genes) > 0) {
        # 检查是否是主成分名
        if (any(grepl("^tt\\.\\d+$", genes) | grepl("^PC\\d+$", genes))) {
          if (verbose) cat("警告：提取到的可能是主成分名而非基因名\n")
          # 对于包含plsRcox的组合，尝试从其他部分提取真实基因名
          if (grepl(" \\+ plsRcox", model_name)) {
            if (verbose) cat("检测到plsRcox组合模型，尝试从主组件提取真实基因名\n")
            # 提取主组件名称
            main_component <- gsub(" \\+ plsRcox", "", model_name)
            if (verbose) cat(sprintf("主组件: %s\n", main_component))
            
            # 根据主组件类型尝试提取
            if (grepl("Lasso", main_component)) {
              main_genes <- try_extract_lasso_component(model_obj, verbose)
            } else if (grepl("CoxBoost", main_component)) {
              main_genes <- try_extract_coxboost_component(model_obj, verbose)
            } else if (grepl("RSF", main_component)) {
              main_genes <- try_extract_rsf_component(model_obj, verbose)
            } else if (grepl("StepCox", main_component)) {
              main_genes <- try_extract_stepcox_component(model_obj, verbose)
            } else {
              main_genes <- NULL
            }
            
            if (!is.null(main_genes) && length(main_genes) > 0) {
              if (verbose) cat(sprintf("从主组件提取到 %d 个真实基因名\n", length(main_genes)))
              return(main_genes)
            }
          }
        }
        return(genes)
      }
    }
    
    # 避免递归：如果是survivalsvm对象，直接按survival-SVM处理
    if (inherits(model_obj, "survivalsvm")) {
      if (verbose) cat("组合模型对象是survivalsvm类型，直接提取基因\n")
      if (!is.null(model_obj$var.names)) {
        if (verbose) cat("从var.names字段提取基因\n")
        return(model_obj$var.names)
      }
    }
    
    genes <- extract_genes_from_ensemble(model_obj, model_name, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      # 检查是否是主成分名
      if (any(grepl("^tt\\.\\d+$", genes) | grepl("^PC\\d+$", genes))) {
        if (verbose) cat("警告：提取到的可能是主成分名而非基因名\n")
        # 对于包含plsRcox的组合，可能需要从其他部分提取基因
        if (grepl(" \\+ plsRcox", model_name)) {
          if (verbose) cat("尝试从组合的其他部分提取真实基因名\n")
          # 这里可以添加更多逻辑，但已经在上面的special_models中处理了
        }
      }
      return(genes)
    }
  }
  
  # 检查是否是SuperPC模型
  if (model_name == "SuperPC") {
    if (verbose) cat("识别为单独的SuperPC模型\n")
    genes <- extract_genes_from_superpc_direct(model_obj, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      return(genes)
    }
  }
  
  # 检查是否是plsRcoxmodel对象（StepCox[both] + plsRcox的情况）
  if (inherits(model_obj, "plsRcoxmodel")) {
    if (verbose) cat("识别为plsRcoxmodel对象\n")
    genes <- extract_genes_from_plsRcoxmodel(model_obj, model_name, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      return(genes)
    }
  }
  
  # 检查是否是SuperPC相关的组合模型
  if (grepl("SuperPC", model_name) && grepl(" \\+ ", model_name)) {
    if (verbose) cat("识别为SuperPC组合模型\n")
    genes <- extract_genes_from_superpc_ensemble(model_obj, model_name, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      return(genes)
    }
  }
  
  # 检查是否是survival-SVM相关的组合模型
  if ((grepl("survival-SVM", model_name) || grepl("survivalsvm", model_name)) && grepl(" \\+ ", model_name)) {
    if (verbose) cat("识别为survival-SVM组合模型\n")
    # 避免递归：如果是survivalsvm对象，直接提取
    if (inherits(model_obj, "survivalsvm")) {
      if (verbose) cat("对象是survivalsvm类型，直接提取基因\n")
      if (!is.null(model_obj$var.names)) {
        if (verbose) cat("从var.names字段提取基因\n")
        return(model_obj$var.names)
      }
    }
    genes <- extract_genes_from_svm_ensemble(model_obj, model_name, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      return(genes)
    }
  }
  
  # 根据模型类型提取基因
  if (grepl("Enet", model_name)) {
    # Elastic Net模型
    if (inherits(model_obj, "cv.glmnet")) {
      if (verbose) cat("识别为: Enet (cv.glmnet)\n")
      coef_mat <- as.matrix(coef(model_obj, s = "lambda.min"))
      genes <- rownames(coef_mat)[which(coef_mat != 0)]
      # 移除可能的截距项
      genes <- genes[!grepl("Intercept", genes)]
      return(genes)
    }
  } else if (grepl("Lasso", model_name) && !grepl("StepCox", model_name) && !grepl("RSF", model_name)) {
    # 单独的Lasso模型
    if (inherits(model_obj, "cv.glmnet")) {
      if (verbose) cat("识别为: Lasso (cv.glmnet)\n")
      coef_mat <- as.matrix(coef(model_obj, s = "lambda.min"))
      genes <- rownames(coef_mat)[which(coef_mat != 0)]
      genes <- genes[!grepl("Intercept", genes)]
      return(genes)
    }
  } else if (grepl("Ridge", model_name) && !grepl("StepCox", model_name) && !grepl("RSF", model_name)) {
    # 单独的Ridge模型
    if (is.list(model_obj) && "fit" %in% names(model_obj)) {
      if (verbose) cat("识别为: Ridge\n")
      coef_mat <- as.matrix(coef(model_obj$fit, s = model_obj$cv.fit$lambda.min))
      genes <- rownames(coef_mat)[which(coef_mat != 0)]
      genes <- genes[!grepl("Intercept", genes)]
      return(genes)
    }
  } else if (grepl("StepCox", model_name)) {
    # StepCox模型
    if (inherits(model_obj, "coxph")) {
      if (verbose) cat("识别为: StepCox (coxph)\n")
      genes <- names(which(!is.na(coef(model_obj)) & coef(model_obj) != 0))
      return(genes)
    } else if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "coxph")) {
      if (verbose) cat("识别为: StepCox (包含fit的列表)\n")
      genes <- names(which(!is.na(coef(model_obj$fit)) & coef(model_obj$fit) != 0))
      return(genes)
    }
  } else if (grepl("CoxBoost", model_name) && !grepl("StepCox", model_name) && !grepl("RSF", model_name) && !grepl("Lasso", model_name)) {
    # CoxBoost模型
    if (inherits(model_obj, "CoxBoost")) {
      if (verbose) cat("识别为: CoxBoost\n")
      coef_mat <- coef(model_obj)
      genes <- names(which(coef_mat != 0))
      return(genes)
    }
  } else if (grepl("plsRcox", model_name)) {
    # plsRcox模型 - 特殊处理
    if (inherits(model_obj, "plsRcox")) {
      if (verbose) cat("识别为: plsRcox\n")
      genes <- extract_plsRcox_genes_safely(model_obj, verbose)
      return(genes)
    } else if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "plsRcox")) {
      if (verbose) cat("识别为: plsRcox (包含fit的列表)\n")
      genes <- extract_plsRcox_genes_safely(model_obj$fit, verbose)
      return(genes)
    }
  } else if (model_name == "RSF" || grepl("^RSF$", model_name)) {
    # RSF模型
    if (inherits(model_obj, "rfsrc")) {
      if (verbose) cat("识别为: RSF (rfsrc)\n")
      var_importance <- model_obj$importance
      return(names(var_importance))
    } else if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "rfsrc")) {
      if (verbose) cat("识别为: RSF (包含fit的列表)\n")
      var_importance <- model_obj$fit$importance
      return(names(var_importance))
    }
  } else if (grepl("GBM", model_name) && !grepl("StepCox", model_name) && !grepl("RSF", model_name) && !grepl("Lasso", model_name)) {
    # GBM模型
    if (is.list(model_obj) && "fit" %in% names(model_obj)) {
      if (verbose) cat("识别为: GBM\n")
      gbm_fit <- model_obj$fit
      rel_inf <- summary(gbm_fit, plotit = FALSE)
      return(as.character(rel_inf$var))
    }
  } else if (grepl("SuperPC", model_name)) {
    # SuperPC模型（单独的或组合中的）
    if (verbose) cat("识别为: SuperPC模型\n")
    genes <- extract_genes_from_superpc_direct(model_obj, verbose)
    if (!is.null(genes) && length(genes) > 0) {
      return(genes)
    }
  } else if (grepl("survival-SVM", model_name) || grepl("survivalsvm", model_name)) {
    # survival-SVM模型
    if (inherits(model_obj, "survivalsvm")) {
      if (verbose) cat("识别为: survival-SVM\n")
      # 直接从var.names字段提取（根据调试结果）
      if (!is.null(model_obj$var.names)) {
        if (verbose) cat("从var.names字段提取基因\n")
        return(model_obj$var.names)
      }
      # 尝试从模型对象中提取特征名
      genes <- try_extract_svm_genes_direct(model_obj, verbose)
      if (!is.null(genes)) {
        return(genes)
      }
      return(NULL)
    } else if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "survivalsvm")) {
      if (verbose) cat("识别为: survival-SVM (包含fit的列表)\n")
      if (!is.null(model_obj$fit$var.names)) {
        if (verbose) cat("从fit$var.names字段提取基因\n")
        return(model_obj$fit$var.names)
      }
      genes <- try_extract_svm_genes_direct(model_obj$fit, verbose)
      if (!is.null(genes)) {
        return(genes)
      }
      return(NULL)
    }
  }
  
  # 对于组合模型，尝试从模型对象中提取
  if (inherits(model_obj, "cv.glmnet")) {
    # 如果是glmnet对象，尝试提取系数
    if (verbose) cat("通用提取: cv.glmnet\n")
    coef_mat <- as.matrix(coef(model_obj, s = "lambda.min"))
    genes <- rownames(coef_mat)[which(coef_mat != 0)]
    genes <- genes[!grepl("Intercept", genes)]
    return(genes)
  } else if (inherits(model_obj, "coxph")) {
    # 如果是coxph对象
    if (verbose) cat("通用提取: coxph\n")
    genes <- names(which(!is.na(coef(model_obj)) & coef(model_obj) != 0))
    return(genes)
  } else if (inherits(model_obj, "CoxBoost")) {
    # 如果是CoxBoost对象
    if (verbose) cat("通用提取: CoxBoost\n")
    coef_mat <- coef(model_obj)
    genes <- names(which(coef_mat != 0))
    return(genes)
  }
  
  # 尝试通用方法提取特征名称
  if (verbose) cat("尝试通用特征提取方法\n")
  genes <- try_extract_features_generic(model_obj, model_name, verbose)
  if (!is.null(genes)) {
    return(genes)
  }
  
  # 对于其他类型，返回NULL
  if (verbose) cat("无法提取基因\n")
  return(NULL)
}

#' 从特殊组合模型中提取基因
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_special_ensemble <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat(sprintf("专门处理特殊组合模型: %s\n", model_name))
  
  # 根据模型名称采取不同的策略
  if (grepl("RSF \\+ SuperPC", model_name)) {
    if (verbose) cat("处理RSF + SuperPC组合\n")
    
    # 检查是否为未命名列表（根据调试结果）
    if (is.list(model_obj) && length(model_obj) == 2 && (is.null(names(model_obj)) || all(names(model_obj) == ""))) {
      if (verbose) cat("模型为未命名列表，尝试按位置提取\n")
      # 第一个元素可能是RSF，第二个可能是SuperPC
      rsf_genes <- try_extract_rsf_component(model_obj[[1]], verbose)
      superpc_genes <- try_extract_superpc_component(model_obj[[2]], verbose)
      
      # 合并基因
      all_genes <- unique(c(rsf_genes, superpc_genes))
      all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
      
      if (length(all_genes) > 0) {
        if (verbose) cat(sprintf("从RSF + SuperPC提取了 %d 个基因\n", length(all_genes)))
        return(all_genes)
      }
    } else {
      # 否则按原来的方法
      if (verbose) cat("使用标准方法提取\n")
      rsf_genes <- try_extract_rsf_component(model_obj, verbose)
      superpc_genes <- try_extract_superpc_component(model_obj, verbose)
      
      # 合并基因，优先使用RSF的基因
      all_genes <- unique(c(rsf_genes, superpc_genes))
      all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
      
      if (length(all_genes) > 0) {
        if (verbose) cat(sprintf("从RSF + SuperPC提取了 %d 个基因\n", length(all_genes)))
        return(all_genes)
      }
    }
    
  } else if (grepl("StepCox\\[", model_name) && grepl(" \\+ RSF", model_name)) {
    if (verbose) cat("处理StepCox + RSF组合\n")
    
    # 提取StepCox类型
    stepcox_type <- gsub(".*(StepCox\\[[^\\]]+\\]).*", "\\1", model_name)
    if (verbose) cat(sprintf("StepCox类型: %s\n", stepcox_type))
    
    # 尝试从StepCox部分提取
    stepcox_genes <- try_extract_stepcox_component(model_obj, verbose)
    
    # 尝试从RSF部分提取
    rsf_genes <- try_extract_rsf_component(model_obj, verbose)
    
    # 合并基因，优先使用StepCox的基因
    all_genes <- unique(c(stepcox_genes, rsf_genes))
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    
    if (length(all_genes) > 0) {
      if (verbose) cat(sprintf("从%s + RSF提取了 %d 个基因\n", stepcox_type, length(all_genes)))
      return(all_genes)
    }
    
  } else if (grepl("Lasso \\+ RSF", model_name)) {
    if (verbose) cat("处理Lasso + RSF组合\n")
    
    # 尝试从Lasso部分提取
    lasso_genes <- try_extract_lasso_component(model_obj, verbose)
    
    # 尝试从RSF部分提取
    rsf_genes <- try_extract_rsf_component(model_obj, verbose)
    
    # 合并基因
    all_genes <- unique(c(lasso_genes, rsf_genes))
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    
    if (length(all_genes) > 0) {
      if (verbose) cat(sprintf("从Lasso + RSF提取了 %d 个基因\n", length(all_genes)))
      return(all_genes)
    }
    
  } else if (grepl("StepCox\\[", model_name) && grepl(" \\+ plsRcox", model_name)) {
    # 专门处理StepCox + plsRcox组合
    if (verbose) cat("处理StepCox + plsRcox组合\n")
    
    # 提取StepCox类型
    stepcox_type <- gsub(".*(StepCox\\[[^\\]]+\\]).*", "\\1", model_name)
    if (verbose) cat(sprintf("StepCox类型: %s\n", stepcox_type))
    
    # 检查模型结构：可能是列表，第一个元素是StepCox，第二个是plsRcox
    if (is.list(model_obj) && length(model_obj) >= 2) {
      if (verbose) cat("模型是列表结构，尝试按位置提取\n")
      
      # 尝试从第一个元素提取StepCox基因
      stepcox_genes <- NULL
      if (!is.null(model_obj[[1]])) {
        stepcox_obj <- model_obj[[1]]
        if (verbose) cat(sprintf("第一个元素类型: %s\n", paste(class(stepcox_obj), collapse=", ")))
        
        # 尝试从StepCox对象提取基因
        if (inherits(stepcox_obj, "coxph")) {
          stepcox_genes <- names(which(!is.na(coef(stepcox_obj)) & coef(stepcox_obj) != 0))
        } else if (is.list(stepcox_obj) && "fit" %in% names(stepcox_obj) && inherits(stepcox_obj$fit, "coxph")) {
          stepcox_genes <- names(which(!is.na(coef(stepcox_obj$fit)) & coef(stepcox_obj$fit) != 0))
        }
      }
      
      # 如果从第一个元素提取失败，尝试通用方法
      if (is.null(stepcox_genes) || length(stepcox_genes) == 0) {
        if (verbose) cat("从第一个元素提取失败，尝试通用StepCox提取\n")
        stepcox_genes <- try_extract_stepcox_component(model_obj, verbose)
      }
      
      if (!is.null(stepcox_genes) && length(stepcox_genes) > 0) {
        if (verbose) cat(sprintf("从StepCox组件提取了 %d 个基因\n", length(stepcox_genes)))
        return(stepcox_genes)
      } else {
        if (verbose) cat("从StepCox部分未提取到基因，尝试从plsRcox部分提取（可能是主成分名）\n")
        plsRcox_genes <- try_extract_plsRcox_component(model_obj, verbose)
        if (!is.null(plsRcox_genes) && length(plsRcox_genes) > 0) {
          if (verbose) cat(sprintf("从plsRcox部分提取了 %d 个特征（可能是主成分名）\n", length(plsRcox_genes)))
          return(plsRcox_genes)
        }
      }
    } else {
      # 如果不是列表结构，使用通用方法
      if (verbose) cat("模型不是列表结构，使用通用方法\n")
      
      # 优先从StepCox部分提取
      stepcox_genes <- try_extract_stepcox_component(model_obj, verbose)
      
      if (!is.null(stepcox_genes) && length(stepcox_genes) > 0) {
        if (verbose) cat(sprintf("从StepCox组件提取了 %d 个基因\n", length(stepcox_genes)))
        return(stepcox_genes)
      } else {
        # 尝试从plsRcox部分提取
        if (verbose) cat("从StepCox部分未提取到基因，尝试从plsRcox部分提取\n")
        plsRcox_genes <- try_extract_plsRcox_component(model_obj, verbose)
        if (!is.null(plsRcox_genes) && length(plsRcox_genes) > 0) {
          if (verbose) cat(sprintf("从plsRcox部分提取了 %d 个特征（可能是主成分名）\n", length(plsRcox_genes)))
          return(plsRcox_genes)
        }
      }
    }
    
  } else if (grepl("Lasso \\+ plsRcox", model_name) || 
             grepl("CoxBoost \\+ plsRcox", model_name) || 
             grepl("RSF \\+ plsRcox", model_name)) {
    
    if (verbose) cat("处理其他包含plsRcox的组合模型\n")
    
    # 确定模型类型
    if (grepl("Lasso \\+ plsRcox", model_name)) {
      main_component <- "Lasso"
      extract_func <- try_extract_lasso_component
    } else if (grepl("CoxBoost \\+ plsRcox", model_name)) {
      main_component <- "CoxBoost"
      extract_func <- try_extract_coxboost_component
    } else if (grepl("RSF \\+ plsRcox", model_name)) {
      main_component <- "RSF"
      extract_func <- try_extract_rsf_component
    } else {
      main_component <- "Unknown"
      extract_func <- NULL
    }
    
    if (!is.null(extract_func)) {
      if (verbose) cat(sprintf("优先从 %s 组件提取真实基因名\n", main_component))
      
      # 检查模型结构：可能是列表，第一个元素是主组件，第二个是plsRcox
      if (is.list(model_obj) && length(model_obj) >= 2) {
        if (verbose) cat("模型是列表结构，尝试从第一个元素提取\n")
        main_genes <- NULL
        
        # 尝试从第一个元素提取
        if (!is.null(model_obj[[1]])) {
          main_obj <- model_obj[[1]]
          # 根据主组件类型提取
          if (main_component == "Lasso" && inherits(main_obj, "cv.glmnet")) {
            coef_mat <- as.matrix(coef(main_obj, s = "lambda.min"))
            main_genes <- rownames(coef_mat)[which(coef_mat != 0)]
            main_genes <- main_genes[!grepl("Intercept", main_genes)]
          } else if (main_component == "CoxBoost" && inherits(main_obj, "CoxBoost")) {
            coef_mat <- coef(main_obj)
            main_genes <- names(which(coef_mat != 0))
          } else if (main_component == "RSF" && inherits(main_obj, "rfsrc")) {
            var_importance <- main_obj$importance
            if (!is.null(var_importance)) {
              main_genes <- names(var_importance)
            }
          }
        }
        
        # 如果从第一个元素提取失败，尝试通用方法
        if (is.null(main_genes) || length(main_genes) == 0) {
          if (verbose) cat("从第一个元素提取失败，尝试通用提取方法\n")
          main_genes <- extract_func(model_obj, verbose)
        }
      } else {
        # 如果不是列表结构，使用通用方法
        main_genes <- extract_func(model_obj, verbose)
      }
      
      if (!is.null(main_genes) && length(main_genes) > 0) {
        if (verbose) cat(sprintf("从%s组件提取了 %d 个基因\n", main_component, length(main_genes)))
        return(main_genes)
      } else {
        # 如果主组件没有提取到基因，尝试从plsRcox部分提取
        if (verbose) cat(sprintf("从%s组件未提取到基因，尝试从plsRcox部分提取\n", main_component))
        plsRcox_genes <- try_extract_plsRcox_component(model_obj, verbose)
        
        if (!is.null(plsRcox_genes) && length(plsRcox_genes) > 0) {
          # 检查是否是主成分名
          if (any(grepl("^tt\\.\\d+$", plsRcox_genes) | grepl("^PC\\d+$", plsRcox_genes))) {
            if (verbose) cat(sprintf("警告：%s + plsRcox提取到的是主成分名，不是原始基因名\n", main_component))
          }
          return(plsRcox_genes)
        }
      }
    }
    
  } else if (grepl("StepCox\\[", model_name) && grepl(" \\+ SuperPC", model_name)) {
    if (verbose) cat("处理StepCox + SuperPC组合\n")
    
    # 提取StepCox类型
    stepcox_type <- gsub(".*(StepCox\\[[^\\]]+\\]).*", "\\1", model_name)
    if (verbose) cat(sprintf("StepCox类型: %s\n", stepcox_type))
    
    # 对于这种结构（fit + cv.fit），fit是superpc对象
    if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "superpc")) {
      if (verbose) cat("模型结构为fit + cv.fit，fit是superpc对象\n")
      # 首先尝试从StepCox部分提取基因
      stepcox_genes <- try_extract_stepcox_component(model_obj, verbose)
      
      # 然后尝试从SuperPC部分提取
      superpc_genes <- try_extract_superpc_component(model_obj, verbose)
      
      # 合并基因，优先使用StepCox的基因
      all_genes <- unique(c(stepcox_genes, superpc_genes))
      all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
      
      if (length(all_genes) > 0) {
        if (verbose) cat(sprintf("从%s + SuperPC提取了 %d 个基因\n", stepcox_type, length(all_genes)))
        return(all_genes)
      }
    } else {
      # 尝试从StepCox部分提取
      stepcox_genes <- try_extract_stepcox_component(model_obj, verbose)
      
      # 尝试从SuperPC部分提取
      superpc_genes <- try_extract_superpc_component(model_obj, verbose)
      
      # 合并基因，优先使用StepCox的基因
      all_genes <- unique(c(stepcox_genes, superpc_genes))
      all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
      
      if (length(all_genes) > 0) {
        if (verbose) cat(sprintf("从%s + SuperPC提取了 %d 个基因\n", stepcox_type, length(all_genes)))
        return(all_genes)
      }
    }
    
  } else if (grepl("CoxBoost \\+ SuperPC", model_name)) {
    if (verbose) cat("处理CoxBoost + SuperPC组合\n")
    
    # 尝试从CoxBoost部分提取
    coxboost_genes <- try_extract_coxboost_component(model_obj, verbose)
    
    # 尝试从SuperPC部分提取
    superpc_genes <- try_extract_superpc_component(model_obj, verbose)
    
    # 合并基因
    all_genes <- unique(c(coxboost_genes, superpc_genes))
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    
    if (length(all_genes) > 0) {
      if (verbose) cat(sprintf("从CoxBoost + SuperPC提取了 %d 个基因\n", length(all_genes)))
      return(all_genes)
    }
    
  } else if (grepl("Lasso \\+ SuperPC", model_name)) {
    if (verbose) cat("处理Lasso + SuperPC组合\n")
    
    # 尝试从Lasso部分提取
    lasso_genes <- try_extract_lasso_component(model_obj, verbose)
    
    # 尝试从SuperPC部分提取
    superpc_genes <- try_extract_superpc_component(model_obj, verbose)
    
    # 合并基因
    all_genes <- unique(c(lasso_genes, superpc_genes))
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    
    if (length(all_genes) > 0) {
      if (verbose) cat(sprintf("从Lasso + SuperPC提取了 %d 个基因\n", length(all_genes)))
      return(all_genes)
    }
  }
  
  if (verbose) cat("无法从特殊组合模型中提取基因\n")
  return(NULL)
}

#' 从SuperPC模型直接提取基因
#'
#' @param model_obj SuperPC模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_superpc_direct <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("直接提取SuperPC模型基因\n")
  
  # 方法1: 检查是否有feature.scores字段（根据调试结果）
  if (is.list(model_obj) && "feature.scores" %in% names(model_obj) && !is.null(model_obj$feature.scores)) {
    feature_scores <- model_obj$feature.scores
    if (!is.null(names(feature_scores))) {
      if (verbose) cat("从feature.scores字段提取基因\n")
      return(names(feature_scores))
    } else if (is.character(feature_scores) && length(feature_scores) > 0) {
      if (verbose) cat("从feature.scores值提取基因\n")
      return(feature_scores)
    }
  }
  
  # 方法2: 检查fit对象中的feature.scores
  if (is.list(model_obj) && "fit" %in% names(model_obj) && is.list(model_obj$fit)) {
    if ("feature.scores" %in% names(model_obj$fit) && !is.null(model_obj$fit$feature.scores)) {
      feature_scores <- model_obj$fit$feature.scores
      if (!is.null(names(feature_scores))) {
        if (verbose) cat("从fit$feature.scores字段提取基因\n")
        return(names(feature_scores))
      } else if (is.character(feature_scores) && length(feature_scores) > 0) {
        if (verbose) cat("从fit$feature.scores值提取基因\n")
        return(feature_scores)
      }
    }
  }
  
  # 方法3: 检查cv.fit对象中的特征
  if (is.list(model_obj) && "cv.fit" %in% names(model_obj) && is.list(model_obj$cv.fit)) {
    cv_fit <- model_obj$cv.fit
    # 检查是否有基因相关的字段
    feature_fields <- c("feature.scores", "features", "featurenames", "gene.names", "var.names")
    for (field in feature_fields) {
      if (field %in% names(cv_fit) && !is.null(cv_fit[[field]])) {
        field_value <- cv_fit[[field]]
        if (!is.null(names(field_value))) {
          if (verbose) cat(sprintf("从cv.fit$%s字段提取基因\n", field))
          return(names(field_value))
        } else if (is.character(field_value) && length(field_value) > 0) {
          if (verbose) cat(sprintf("从cv.fit$%s值提取基因\n", field))
          return(field_value)
        }
      }
    }
  }
  
  # 方法4: 检查是否有featurenames字段
  if (is.list(model_obj) && "featurenames" %in% names(model_obj) && !is.null(model_obj$featurenames)) {
    if (verbose) cat("从featurenames字段提取基因\n")
    return(model_obj$featurenames)
  }
  
  # 方法5: 检查fit对象中的featurenames
  if (is.list(model_obj) && "fit" %in% names(model_obj) && is.list(model_obj$fit)) {
    if ("featurenames" %in% names(model_obj$fit) && !is.null(model_obj$fit$featurenames)) {
      if (verbose) cat("从fit$featurenames字段提取基因\n")
      return(model_obj$fit$featurenames)
    }
  }
  
  # 方法6: 检查loadings矩阵
  if (is.list(model_obj) && "loadings" %in% names(model_obj) && !is.null(model_obj$loadings)) {
    loadings_mat <- model_obj$loadings
    if (is.matrix(loadings_mat) && !is.null(rownames(loadings_mat))) {
      if (verbose) cat("从loadings矩阵提取基因\n")
      return(rownames(loadings_mat))
    }
  }
  
  # 方法7: 检查rotation矩阵
  if (is.list(model_obj) && "rotation" %in% names(model_obj) && !is.null(model_obj$rotation)) {
    rotation_mat <- model_obj$rotation
    if (is.matrix(rotation_mat) && !is.null(rownames(rotation_mat))) {
      if (verbose) cat("从rotation矩阵提取基因\n")
      return(rownames(rotation_mat))
    }
  }
  
  # 方法8: 检查v矩阵（主成分）
  if (is.list(model_obj) && "v" %in% names(model_obj) && !is.null(model_obj$v)) {
    v_mat <- model_obj$v
    if (is.matrix(v_mat) && !is.null(rownames(v_mat))) {
      if (verbose) cat("从v矩阵提取基因\n")
      return(rownames(v_mat))
    }
  }
  
  # 方法9: 检查训练数据
  if (is.list(model_obj)) {
    # 检查是否有训练数据字段
    train_fields <- c("x", "X", "data", "dataX", "train.x", "trainX")
    for (field in train_fields) {
      if (field %in% names(model_obj) && !is.null(model_obj[[field]])) {
        data_obj <- model_obj[[field]]
        
        if (is.matrix(data_obj) && !is.null(colnames(data_obj))) {
          if (verbose) cat(sprintf("从训练数据字段 %s 提取基因\n", field))
          return(colnames(data_obj))
        } else if (is.data.frame(data_obj) && length(data_obj) > 0) {
          if (verbose) cat(sprintf("从训练数据字段 %s 提取基因\n", field))
          return(names(data_obj))
        }
      }
    }
  }
  
  # 方法10: 对于SuperPC模型的特殊结构（fit为superpc对象）
  if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "superpc")) {
    if (verbose) cat("SuperPC模型结构，尝试从原始特征选择中获取基因名\n")
    # 尝试从原始数据中获取基因名 - 这可能需要外部信息
    if (!is.null(model_obj$call)) {
      call_str <- deparse(model_obj$call)
      if (grepl("featurenames", call_str)) {
        # 尝试从调用中解析featurenames参数
        if (verbose) cat("调用中包含featurenames参数，尝试解析\n")
        # 这里可以进一步解析，但比较复杂
      }
    }
  }
  
  if (verbose) cat("无法从SuperPC模型直接提取基因\n")
  return(NULL)
}

#' 安全地从plsRcox对象提取基因
#'
#' @param model_obj plsRcox对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_plsRcox_genes_safely <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("安全提取plsRcox基因\n")
  
  # 优先从dataX字段提取原始基因名（根据调试信息）
  if (!is.null(model_obj$dataX)) {
    if (is.data.frame(model_obj$dataX)) {
      if (verbose) cat("从dataX字段提取基因名\n")
      return(names(model_obj$dataX))
    } else if (is.matrix(model_obj$dataX) && !is.null(colnames(model_obj$dataX))) {
      if (verbose) cat("从dataX矩阵提取基因名\n")
      return(colnames(model_obj$dataX))
    }
  }
  
  # 尝试提取变量名
  if (!is.null(model_obj$var.names)) {
    if (verbose) cat("从var.names字段提取\n")
    return(model_obj$var.names)
  }
  
  # 尝试提取Xnames
  if (!is.null(model_obj$Xnames)) {
    if (verbose) cat("从Xnames字段提取\n")
    return(model_obj$Xnames)
  }
  
  # 尝试提取系数矩阵的行名
  if (!is.null(model_obj$Coeffs)) {
    coefs <- model_obj$Coeffs
    if (is.matrix(coefs) && !is.null(rownames(coefs))) {
      if (verbose) cat("从Coeffs矩阵提取\n")
      return(rownames(coefs))
    }
  }
  
  # 尝试提取tt矩阵的行名
  if (!is.null(model_obj$tt)) {
    tt_mat <- model_obj$tt
    if (is.matrix(tt_mat) && !is.null(rownames(tt_mat))) {
      if (verbose) cat("从tt矩阵提取\n")
      return(rownames(tt_mat))
    }
  }
  
  # 检查是否有loadings矩阵
  if (!is.null(model_obj$loadings)) {
    loadings_mat <- model_obj$loadings
    if (is.matrix(loadings_mat) && !is.null(rownames(loadings_mat))) {
      if (verbose) cat("从loadings矩阵提取\n")
      return(rownames(loadings_mat))
    }
  }
  
  # 检查ExpliX字段（解释变量矩阵）
  if (!is.null(model_obj$ExpliX)) {
    expliX <- model_obj$ExpliX
    if (is.matrix(expliX) && !is.null(colnames(expliX))) {
      if (verbose) cat("从ExpliX矩阵提取基因名\n")
      return(colnames(expliX))
    }
  }
  
  # 检查XXwotNA字段（去除NA的矩阵）
  if (!is.null(model_obj$XXwotNA)) {
    xx <- model_obj$XXwotNA
    if (is.matrix(xx) && !is.null(colnames(xx))) {
      if (verbose) cat("从XXwotNA矩阵提取基因名\n")
      return(colnames(xx))
    }
  }
  
  return(NULL)
}


# 辅助函数定义
#' 在对象中查找指定类的子对象
#'
#' @param obj 对象
#' @param class_name 类名
#' @param max_depth 最大递归深度
#' @return 找到的对象或NULL
find_object_by_class <- function(obj, class_name, max_depth = 3) {
  if (max_depth < 0) return(NULL)
  
  # 如果对象本身就是指定类
  if (inherits(obj, class_name)) {
    return(obj)
  }
  
  # 如果是列表，递归查找
  if (is.list(obj) && length(obj) > 0) {
    for (i in seq_along(obj)) {
      elem <- obj[[i]]
      result <- find_object_by_class(elem, class_name, max_depth - 1)
      if (!is.null(result)) {
        return(result)
      }
    }
  }
  
  return(NULL)
}

#' 尝试从RSF组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_rsf_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取RSF组件基因\n")
  
  # 如果是列表且长度为2且未命名（针对RSF + SuperPC的特殊情况）
  if (is.list(model_obj) && length(model_obj) == 2 && (is.null(names(model_obj)) || all(names(model_obj) == ""))) {
    if (verbose) cat("对象是未命名列表，尝试提取第一个元素作为RSF\n")
    # 尝试第一个元素
    model_obj <- model_obj[[1]]
  }
  
  # 如果是列表，检查是否是包含RSF的列表结构
  if (is.list(model_obj)) {
    # 检查是否是列表的第一个元素
    if (length(model_obj) >= 1 && !is.null(model_obj[[1]]) && inherits(model_obj[[1]], "rfsrc")) {
      if (verbose) cat("从第一个元素找到rfsrc对象\n")
      var_importance <- model_obj[[1]]$importance
      if (!is.null(var_importance)) {
        return(names(var_importance))
      }
    }
    
    # 查找rfsrc对象
    rfsrc_obj <- find_object_by_class(model_obj, "rfsrc")
    if (!is.null(rfsrc_obj)) {
      if (verbose) cat("找到RSF对象\n")
      var_importance <- rfsrc_obj$importance
      if (!is.null(var_importance)) {
        return(names(var_importance))
      }
    }
    
    # 查找importance字段
    if ("importance" %in% names(model_obj) && !is.null(model_obj$importance)) {
      if (verbose) cat("从importance字段提取RSF基因\n")
      return(names(model_obj$importance))
    }
  }
  
  return(NULL)
}

#' 尝试从StepCox组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_stepcox_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取StepCox组件基因\n")
  
  # 检查模型结构：如果是列表，可能是组合模型
  if (is.list(model_obj)) {
    # 检查是否是包含StepCox的列表结构
    if (length(model_obj) >= 2 && !is.null(names(model_obj))) {
      # 检查是否有明显的StepCox相关字段
      if ("fit" %in% names(model_obj) && inherits(model_obj$fit, "coxph")) {
        if (verbose) cat("从fit字段找到coxph对象\n")
        genes <- names(which(!is.na(coef(model_obj$fit)) & coef(model_obj$fit) != 0))
        if (length(genes) > 0) {
          return(genes)
        }
      }
    }
    
    # 尝试查找coxph对象
    coxph_obj <- find_object_by_class(model_obj, "coxph")
    if (!is.null(coxph_obj)) {
      if (verbose) cat("找到StepCox对象\n")
      genes <- names(which(!is.na(coef(coxph_obj)) & coef(coxph_obj) != 0))
      if (length(genes) > 0) {
        return(genes)
      }
    }
    
    # 检查是否有系数字段
    if ("coefficients" %in% names(model_obj) && !is.null(model_obj$coefficients)) {
      coefs <- model_obj$coefficients
      if (!is.null(names(coefs))) {
        genes <- names(which(coefs != 0))
        if (length(genes) > 0) {
          if (verbose) cat("从coefficients字段提取StepCox基因\n")
          return(genes)
        }
      }
    }
    
    # 检查是否是列表的第一个元素
    if (length(model_obj) >= 1 && !is.null(model_obj[[1]])) {
      if (verbose) cat("尝试从第一个元素提取\n")
      return(try_extract_stepcox_component(model_obj[[1]], verbose))
    }
  }
  
  return(NULL)
}

#' 尝试从CoxBoost组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_coxboost_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取CoxBoost组件基因\n")
  
  # 查找CoxBoost对象
  coxboost_obj <- find_object_by_class(model_obj, "CoxBoost")
  if (!is.null(coxboost_obj)) {
    if (verbose) cat("找到CoxBoost对象\n")
    coef_mat <- coef(coxboost_obj)
    genes <- names(which(coef_mat != 0))
    if (length(genes) > 0) {
      return(genes)
    }
  }
  
  return(NULL)
}

#' 尝试从Lasso组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_lasso_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取Lasso组件基因\n")
  
  # 查找cv.glmnet对象
  glmnet_obj <- find_object_by_class(model_obj, "cv.glmnet")
  if (!is.null(glmnet_obj)) {
    if (verbose) cat("找到Lasso对象\n")
    coef_mat <- as.matrix(coef(glmnet_obj, s = "lambda.min"))
    genes <- rownames(coef_mat)[which(coef_mat != 0)]
    genes <- genes[!grepl("Intercept", genes)]
    if (length(genes) > 0) {
      return(genes)
    }
  }
  
  return(NULL)
}

#' 尝试从SVM组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_svm_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取SVM组件基因\n")
  
  # 查找survivalsvm对象
  svm_obj <- find_object_by_class(model_obj, "survivalsvm")
  if (!is.null(svm_obj)) {
    if (verbose) cat("找到survivalsvm对象\n")
    if (!is.null(svm_obj$var.names)) {
      if (verbose) cat("从var.names字段提取SVM基因\n")
      return(svm_obj$var.names)
    }
    return(try_extract_svm_genes_direct(svm_obj, verbose))
  }
  
  # 尝试从训练数据中提取
  if (is.list(model_obj)) {
    train_fields <- c("x", "X", "data", "train.x", "trainX")
    for (field in train_fields) {
      if (field %in% names(model_obj) && !is.null(model_obj[[field]])) {
        data_obj <- model_obj[[field]]
        if (is.matrix(data_obj) && !is.null(colnames(data_obj))) {
          if (verbose) cat(sprintf("从训练数据字段 %s 提取SVM基因\n", field))
          return(colnames(data_obj))
        }
      }
    }
  }
  
  return(NULL)
}

#' 尝试从plsRcox组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_plsRcox_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取plsRcox组件基因\n")
  
  # 查找plsRcox对象
  plsRcox_obj <- find_object_by_class(model_obj, "plsRcox")
  if (!is.null(plsRcox_obj)) {
    if (verbose) cat("找到plsRcox对象\n")
    genes <- extract_plsRcox_genes_safely(plsRcox_obj, verbose)
    return(genes)
  }
  
  # 查找plsRcoxmodel对象
  plsRcoxmodel_obj <- find_object_by_class(model_obj, "plsRcoxmodel")
  if (!is.null(plsRcoxmodel_obj)) {
    if (verbose) cat("找到plsRcoxmodel对象\n")
    genes <- extract_genes_from_plsRcoxmodel(plsRcoxmodel_obj, "plsRcoxmodel", verbose)
    return(genes)
  }
  
  return(NULL)
}

#' 从plsRcoxmodel对象提取基因
#'
#' @param model_obj plsRcoxmodel对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_plsRcoxmodel <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat("从plsRcoxmodel提取基因\n")
  
  # 首先检查dataX字段（包含原始基因名）
  if (is.list(model_obj) && "dataX" %in% names(model_obj)) {
    dataX <- model_obj$dataX
    if (is.data.frame(dataX)) {
      genes <- names(dataX)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从dataX字段提取了 %d 个基因\n", length(genes)))
        return(genes)
      }
    } else if (is.matrix(dataX) && !is.null(colnames(dataX))) {
      genes <- colnames(dataX)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从dataX矩阵提取了 %d 个基因\n", length(genes)))
        return(genes)
      }
    }
  }
  
  # 检查ExpliX字段
  if (is.list(model_obj) && "ExpliX" %in% names(model_obj)) {
    expliX <- model_obj$ExpliX
    if (is.matrix(expliX) && !is.null(colnames(expliX))) {
      genes <- colnames(expliX)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从ExpliX矩阵提取了 %d 个基因\n", length(genes)))
        return(genes)
      }
    }
  }
  
  # 检查XXwotNA字段
  if (is.list(model_obj) && "XXwotNA" %in% names(model_obj)) {
    xx <- model_obj$XXwotNA
    if (is.matrix(xx) && !is.null(colnames(xx))) {
      genes <- colnames(xx)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从XXwotNA矩阵提取了 %d 个基因\n", length(genes)))
        return(genes)
      }
    }
  }
  
  # 检查系数矩阵
  if (is.list(model_obj) && "Coeffs" %in% names(model_obj)) {
    coefs <- model_obj$Coeffs
    if (is.matrix(coefs) && !is.null(rownames(coefs))) {
      genes <- rownames(coefs)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从Coeffs矩阵提取了 %d 个基因\n", length(genes)))
        return(genes)
      }
    }
  }
  
  # 检查变量名
  if (is.list(model_obj) && "var.names" %in% names(model_obj)) {
    var_names <- model_obj$var.names
    if (is.character(var_names) && length(var_names) > 0) {
      if (verbose) cat(sprintf("从var.names字段提取了 %d 个基因\n", length(var_names)))
      return(var_names)
    }
  }
  
  # 检查tt矩阵的行名
  if (is.list(model_obj) && "tt" %in% names(model_obj)) {
    tt_mat <- model_obj$tt
    if (is.matrix(tt_mat) && !is.null(rownames(tt_mat))) {
      genes <- rownames(tt_mat)
      if (length(genes) > 0) {
        if (verbose) cat(sprintf("从tt矩阵提取了 %d 个基因\n", length(genes)))
        # 检查是否是主成分名
        if (any(grepl("^tt\\.\\d+$", genes) | grepl("^PC\\d+$", genes))) {
          if (verbose) cat("警告：提取到的是主成分名，不是原始基因名\n")
        }
        return(genes)
      }
    }
  }
  
  # 检查FinalModel字段中的原始基因信息
  if (is.list(model_obj) && "FinalModel" %in% names(model_obj)) {
    final_model <- model_obj$FinalModel
    if (inherits(final_model, "coxph")) {
      # 检查是否有原始数据的引用
      if (!is.null(final_model$model)) {
        model_data <- final_model$model
        # 尝试从模型数据中提取变量名
        if (is.data.frame(model_data) && length(model_data) > 0) {
          # 排除时间和状态列
          potential_genes <- names(model_data)
          # 移除可能的生存时间/状态列
          time_status_cols <- c("time", "status", "surv.time", "surv.status", 
                               "Time", "Status", "OS.time", "OS")
          genes <- setdiff(potential_genes, time_status_cols)
          if (length(genes) > 0) {
            if (verbose) cat(sprintf("从FinalModel模型数据提取了 %d 个基因\n", length(genes)))
            return(genes)
          }
        }
      }
    }
  }
  
  return(NULL)
}


#' 尝试从SuperPC组件提取基因
#'
#' @param model_obj 模型对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_superpc_component <- function(model_obj, verbose = FALSE) {
  if (verbose) cat("尝试提取SuperPC组件基因\n")
  
  # 如果是列表且长度为2且未命名（针对RSF + SuperPC的特殊情况）
  if (is.list(model_obj) && length(model_obj) == 2 && (is.null(names(model_obj)) || all(names(model_obj) == ""))) {
    if (verbose) cat("对象是未命名列表，尝试提取第二个元素作为SuperPC\n")
    # 尝试第二个元素
    model_obj <- model_obj[[2]]
  }
  
  # 使用直接提取SuperPC基因的方法
  genes <- extract_genes_from_superpc_direct(model_obj, verbose)
  
  # 如果没找到基因，检查是否是包含fit和cv.fit的列表结构
  if (is.null(genes) || length(genes) == 0) {
    if (is.list(model_obj) && "fit" %in% names(model_obj) && inherits(model_obj$fit, "superpc")) {
      if (verbose) cat("模型结构为fit + cv.fit，尝试从fit对象提取\n")
      genes <- extract_genes_from_superpc_direct(model_obj$fit, verbose)
    }
  }
  
  return(genes)
}

#' 尝试直接从SVM对象提取基因
#'
#' @param svm_obj survivalsvm对象
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_svm_genes_direct <- function(svm_obj, verbose = FALSE) {
  if (verbose) cat("尝试直接提取SVM基因\n")
  
  # 方法1: 检查模型拟合对象
  if (is.list(svm_obj) && "model.fit" %in% names(svm_obj)) {
    model_fit <- svm_obj$model.fit
    if (is.list(model_fit)) {
      # 检查是否有特征权重
      if ("weights" %in% names(model_fit) && !is.null(model_fit$weights)) {
        weights <- model_fit$weights
        if (!is.null(names(weights))) {
          if (verbose) cat("从model.fit$weights提取基因\n")
          return(names(weights))
        }
      }
    }
  }
  
  # 方法2: 检查训练数据
  if (is.list(svm_obj) && "x" %in% names(svm_obj)) {
    x_data <- svm_obj$x
    if (is.matrix(x_data) && !is.null(colnames(x_data))) {
      if (verbose) cat("从x字段提取基因\n")
      return(colnames(x_data))
    }
  }
  
  # 方法3: 检查调用参数
  if (!is.null(svm_obj$call)) {
    call_str <- deparse(svm_obj$call)
    # 尝试解析调用中的x参数
    if (grepl("x\\s*=", call_str)) {
      if (verbose) cat("尝试从调用中解析基因名\n")
      # 这里可以进一步解析，但比较复杂
    }
  }
  
  if (verbose) cat("无法直接提取SVM基因\n")
  return(NULL)
}

#' 从组合模型中提取基因
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_ensemble <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat("从组合模型中提取基因\n")
  
  # 检查是否是包含多个模型的列表
  if (is.list(model_obj)) {
    # 首先尝试直接提取基因（对于已经是完整模型对象的情况）
    if (inherits(model_obj, "survivalsvm") && !is.null(model_obj$var.names)) {
      if (verbose) cat("组合模型对象是survivalsvm类型，直接从var.names提取\n")
      return(model_obj$var.names)
    }
    
    # 尝试从每个组件中提取基因
    all_genes <- c()
    processed_count <- 0
    
    for (i in seq_along(model_obj)) {
      elem_name <- names(model_obj)[i]
      if (is.null(elem_name) || elem_name == "") {
        elem_name <- paste0("[[", i, "]]")
      }
      
      elem <- model_obj[[i]]
      
      if (verbose) {
        cat(sprintf("  处理元素 %s: ", elem_name))
        cat(sprintf("类型: %s\n", paste(class(elem), collapse=", ")))
      }
      
      # 跳过NULL元素
      if (is.null(elem)) {
        if (verbose) cat("    -> 元素为NULL，跳过\n")
        next
      }
      
      # 检查元素是否与父对象相同（避免递归调用自身）
      if (identical(elem, model_obj)) {
        if (verbose) cat("    -> 元素与父对象相同，跳过\n")
        next
      }
      
      # 对于已知的模型类型，使用简化提取方法
      if (inherits(elem, "survivalsvm")) {
        if (!is.null(elem$var.names)) {
          if (verbose) cat(sprintf("    -> 从survivalsvm的var.names提取 %d 个基因\n", length(elem$var.names)))
          all_genes <- c(all_genes, elem$var.names)
          processed_count <- processed_count + 1
          next
        }
      } else if (inherits(elem, "cv.glmnet")) {
        tryCatch({
          coef_mat <- as.matrix(coef(elem, s = "lambda.min"))
          genes <- rownames(coef_mat)[which(coef_mat != 0)]
          genes <- genes[!grepl("Intercept", genes)]
          if (length(genes) > 0) {
            if (verbose) cat(sprintf("    -> 从cv.glmnet提取 %d 个基因\n", length(genes)))
            all_genes <- c(all_genes, genes)
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取cv.glmnet基因时出错: %s\n", e$message))
        })
        next
      } else if (inherits(elem, "coxph")) {
        tryCatch({
          genes <- names(which(!is.na(coef(elem)) & coef(elem) != 0))
          if (length(genes) > 0) {
            if (verbose) cat(sprintf("    -> 从coxph提取 %d 个基因\n", length(genes)))
            all_genes <- c(all_genes, genes)
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取coxph基因时出错: %s\n", e$message))
        })
        next
      } else if (inherits(elem, "CoxBoost")) {
        tryCatch({
          coef_mat <- coef(elem)
          genes <- names(which(coef_mat != 0))
          if (length(genes) > 0) {
            if (verbose) cat(sprintf("    -> 从CoxBoost提取 %d 个基因\n", length(genes)))
            all_genes <- c(all_genes, genes)
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取CoxBoost基因时出错: %s\n", e$message))
        })
        next
      } else if (inherits(elem, "rfsrc")) {
        tryCatch({
          var_importance <- elem$importance
          if (!is.null(var_importance) && !is.null(names(var_importance))) {
            if (verbose) cat(sprintf("    -> 从rfsrc提取 %d 个基因\n", length(names(var_importance))))
            all_genes <- c(all_genes, names(var_importance))
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取rfsrc基因时出错: %s\n", e$message))
        })
        next
      } else if (inherits(elem, "superpc") || inherits(elem, "superpc.cv")) {
        # 对于superpc和superpc.cv对象，使用改进的提取方法
        genes <- extract_genes_from_superpc_direct(elem, FALSE)
        if (!is.null(genes) && length(genes) > 0) {
          if (verbose) cat(sprintf("    -> 从superpc提取 %d 个基因\n", length(genes)))
          all_genes <- c(all_genes, genes)
          processed_count <- processed_count + 1
        } else {
          if (verbose) cat("    -> 从superpc对象未提取到基因\n")
        }
        next
      } else if (inherits(elem, "gbm")) {
        tryCatch({
          rel_inf <- summary(elem, plotit = FALSE)
          if (!is.null(rel_inf$var)) {
            genes <- as.character(rel_inf$var)
            if (length(genes) > 0) {
              if (verbose) cat(sprintf("    -> 从gbm提取 %d 个基因\n", length(genes)))
              all_genes <- c(all_genes, genes)
              processed_count <- processed_count + 1
            }
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取gbm基因时出错: %s\n", e$message))
        })
        next
      } else if (inherits(elem, "plsRcox") || inherits(elem, "plsRcoxmodel")) {
        tryCatch({
          genes <- extract_plsRcox_genes_safely(elem, FALSE)
          if (!is.null(genes) && length(genes) > 0) {
            if (verbose) cat(sprintf("    -> 从plsRcox提取 %d 个基因\n", length(genes)))
            all_genes <- c(all_genes, genes)
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 提取plsRcox基因时出错: %s\n", e$message))
        })
        next
      } else if (is.list(elem) && length(elem) > 0) {
        # 对于其他列表类型元素，尝试递归提取
        if (verbose) cat("    -> 尝试递归提取列表元素\n")
        tryCatch({
          # 使用深度限制避免无限递归
          sub_genes <- extract_genes_from_model_recursive(elem, model_name, verbose, max_depth = 1)
          if (!is.null(sub_genes) && length(sub_genes) > 0) {
            if (verbose) cat(sprintf("    -> 递归提取到 %d 个基因\n", length(sub_genes)))
            all_genes <- c(all_genes, sub_genes)
            processed_count <- processed_count + 1
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("    -> 递归提取时出错: %s\n", e$message))
        })
        next
      }
      
      if (verbose) cat("    -> 无法处理此元素类型\n")
    }
    
    if (verbose) cat(sprintf("总计处理了 %d 个元素，提取到 %d 个唯一基因\n", processed_count, length(unique(all_genes))))
    
    if (length(all_genes) > 0) {
      all_genes <- unique(all_genes)
      return(all_genes)
    }
  }
  
  return(NULL)
}

#' 递归提取基因（带深度限制）
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @param max_depth 最大递归深度
#' @param current_depth 当前深度
#' @return 基因名称向量
extract_genes_from_model_recursive <- function(model_obj, model_name, verbose = FALSE, max_depth = 2, current_depth = 0) {
  if (current_depth >= max_depth) {
    if (verbose) cat(sprintf("达到最大递归深度 %d，停止递归\n", max_depth))
    return(NULL)
  }
  
  # 简化版提取，避免无限递归
  if (inherits(model_obj, "survivalsvm") && !is.null(model_obj$var.names)) {
    return(model_obj$var.names)
  } else if (inherits(model_obj, "cv.glmnet")) {
    tryCatch({
      coef_mat <- as.matrix(coef(model_obj, s = "lambda.min"))
      genes <- rownames(coef_mat)[which(coef_mat != 0)]
      genes <- genes[!grepl("Intercept", genes)]
      return(genes)
    }, error = function(e) return(NULL))
  } else if (inherits(model_obj, "coxph")) {
    tryCatch({
      genes <- names(which(!is.na(coef(model_obj)) & coef(model_obj) != 0))
      return(genes)
    }, error = function(e) return(NULL))
  } else if (inherits(model_obj, "CoxBoost")) {
    tryCatch({
      coef_mat <- coef(model_obj)
      genes <- names(which(coef_mat != 0))
      return(genes)
    }, error = function(e) return(NULL))
  } else if (inherits(model_obj, "rfsrc")) {
    tryCatch({
      var_importance <- model_obj$importance
      if (!is.null(var_importance)) {
        return(names(var_importance))
      }
    }, error = function(e) return(NULL))
  } else if (inherits(model_obj, "superpc") || inherits(model_obj, "superpc.cv")) {
    return(extract_genes_from_superpc_direct(model_obj, FALSE))
  } else if (inherits(model_obj, "gbm")) {
    tryCatch({
      rel_inf <- summary(model_obj, plotit = FALSE)
      return(as.character(rel_inf$var))
    }, error = function(e) return(NULL))
  } else if (inherits(model_obj, "plsRcox") || inherits(model_obj, "plsRcoxmodel")) {
    return(extract_plsRcox_genes_safely(model_obj, FALSE))
  } else if (is.list(model_obj)) {
    # 对于列表，尝试提取常见的基因字段
    feature_fields <- c("features", "featurenames", "featureNames", "geneNames", "var.names", 
                        "feature.scores", "Xnames")
    for (field in feature_fields) {
      if (field %in% names(model_obj)) {
        field_value <- model_obj[[field]]
        if (is.character(field_value) && length(field_value) > 0) {
          return(field_value)
        } else if (!is.null(names(field_value))) {
          return(names(field_value))
        }
      }
    }
    
    # 检查训练数据字段
    data_fields <- c("x", "X", "data", "train.x", "trainX")
    for (field in data_fields) {
      if (field %in% names(model_obj)) {
        data_obj <- model_obj[[field]]
        if (is.matrix(data_obj) && !is.null(colnames(data_obj))) {
          return(colnames(data_obj))
        } else if (is.data.frame(data_obj) && length(data_obj) > 0) {
          return(names(data_obj))
        }
      }
    }
  }
  
  return(NULL)
}

#' 从SuperPC组合模型中提取基因
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_superpc_ensemble <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat("从SuperPC组合模型中提取基因\n")
  
  # 首先尝试直接提取
  genes <- extract_genes_from_superpc_direct(model_obj, verbose)
  if (!is.null(genes) && length(genes) > 0) {
    return(genes)
  }
  
  # 尝试从组合的其他组件中提取
  return(extract_genes_from_ensemble(model_obj, model_name, verbose))
}

#' 从SVM组合模型中提取基因
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
extract_genes_from_svm_ensemble <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat("从SVM组合模型中提取基因\n")
  
  # 首先尝试直接提取SVM基因
  if (inherits(model_obj, "survivalsvm")) {
    return(extract_genes_from_model(model_obj, model_name, verbose))
  }
  
  # 尝试从组合的其他组件中提取
  return(extract_genes_from_ensemble(model_obj, model_name, verbose))
}

#' 尝试通用方法提取特征
#'
#' @param model_obj 模型对象
#' @param model_name 模型名称
#' @param verbose 是否显示调试信息
#' @return 基因名称向量
try_extract_features_generic <- function(model_obj, model_name, verbose = FALSE) {
  if (verbose) cat("尝试通用特征提取方法\n")
  
  # 检查是否有变量名字段
  feature_fields <- c("features", "featurenames", "featureNames", "geneNames", "var.names", 
                      "variable.names", "predictorNames", "colnames", "names", "Xnames", "feature.scores")
  
  for (field in feature_fields) {
    if (is.list(model_obj) && field %in% names(model_obj)) {
      field_value <- model_obj[[field]]
      if (is.character(field_value) && length(field_value) > 0) {
        if (verbose) cat(sprintf("从 %s 字段提取基因\n", field))
        return(field_value)
      } else if (!is.null(names(field_value))) {
        if (verbose) cat(sprintf("从 %s 字段的名称提取基因\n", field))
        return(names(field_value))
      }
    }
  }
  
  # 检查是否有训练数据
  data_fields <- c("x", "X", "data", "train.x", "trainX")
  for (field in data_fields) {
    if (is.list(model_obj) && field %in% names(model_obj)) {
      data_obj <- model_obj[[field]]
      if (is.matrix(data_obj) && !is.null(colnames(data_obj))) {
        if (verbose) cat(sprintf("从训练数据字段 %s 提取基因\n", field))
        return(colnames(data_obj))
      } else if (is.data.frame(data_obj) && length(data_obj) > 0) {
        if (verbose) cat(sprintf("从训练数据字段 %s 提取基因\n", field))
        return(names(data_obj))
      }
    }
  }
  
  # 检查是否有系数
  if (is.list(model_obj) && "coefficients" %in% names(model_obj)) {
    coefs <- model_obj$coefficients
    if (!is.null(names(coefs))) {
      genes <- names(coefs)[which(coefs != 0)]
      genes <- genes[!grepl("Intercept", genes)]
      if (length(genes) > 0) {
        if (verbose) cat("从coefficients字段提取基因\n")
        return(genes)
      }
    }
  }
  
  return(NULL)
}

# 使用示例：
# all_genes <- get_genes_by_model(res, verbose = TRUE, debug_missing = TRUE)
# save_genes_to_excel(all_genes, "all_model_genes.xlsx")
