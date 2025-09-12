library(ggplot2)
library(zoo)
library(evd)
library(gridExtra)
library(tibble)

load("run1.RData")
load("run2.RData")
load("run3.RData")
load("run4.RData")


## mean level trend check ##
df <- data.frame(
  time = 1:length(run1[1,1,]),
  value = run1[1,1,]
)

ggplot(df, aes(x = time, y = value)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.1, color = "blue", se = FALSE) +
  labs(title = "Time Series with LOESS Smoothing",
       x = "Time",
       y = "Value") +
  theme_minimal()




## extreme level trend check ##
ts_data <- run1[1,1,]  

# long term trend
window_size <- 365*5  # e.g., 5-year window

# Rolling quantile computation
rolling_q_90 <- rollapply(ts_data, width = window_size,
                          FUN = function(x) quantile(x, probs = 0.9),
                          by = 1, align = "right", fill = NA)

rolling_q_99 <- rollapply(ts_data, width = window_size,
                          FUN = function(x) quantile(x, probs = 0.99),
                          by = 1, align = "right", fill = NA)

rolling_q_999 <- rollapply(ts_data, width = window_size,
                           FUN = function(x) quantile(x, probs = 0.999),
                           by = 1, align = "right", fill = NA)

# Time index for plotting
time_index <- seq_along(rolling_q_90)

# Plot
par(mar=c(5,4.5,4,2)+0.1)
plot(time_index, ts_data, pch = 16, col = "gray", main = "", 
     xlab="Time (date)", ylab="Precipitation (Leadbetters)", cex.lab=1.2)
lines(time_index, rolling_q_90, type = 'l', col = "red", lwd=2)
lines(time_index, rolling_q_99, type = 'l', col = "blue", lwd=2)
lines(time_index, rolling_q_999, type = 'l', col = "green", lwd=2)

legend("topright", 
       legend = c(expression(tau == 0.999), 
                  expression(tau == 0.99), 
                  expression(tau == 0.9)),
       col = c("green", "blue", "red"), 
       lwd = 2, cex=1.3)


## seasonality check ##

window_size <- 180  # e.g., 180-day window

rolling_q_99 <- rollapply(ts_data, width = window_size,
                          FUN = function(x) quantile(x, probs = 0.99),
                          by = 1, align = "right", fill = NA)

period <- c((365*150+1) : (365*160)) # data from 2000 to 2009

plot(time_index[period], ts_data[period], pch = 16, col = "gray", main = "", 
     xlab="Time (date)", ylab="Precipitation (Leadbetters)", cex.lab=1.2)

lines(time_index[period], rolling_q_99[period], type = 'l', col = "blue", lwd=2)


ndays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
abline(v=365*c(150:159)+sum(ndays[1:4]), lty=1, col="red")
abline(v=365*c(150:160)+sum(ndays[1:10]), lty=2)

legend("topright",
       legend = c(expression(tau == 0.99)),
       col = c("blue"),
       lwd = 2, cex=1.3)



ndays.cumsum <- cumsum(ndays)
month.index <- sapply(1:365, \(i) which.max(i <= ndays.cumsum))

## generate 60225 x 25 matrix
## the first row: (1,1), (2,1), (3,1), (4,1), (5,1), (2,1) ...
collapse_grid <- function(RUN){
  X <- apply(RUN,3,c) |> t()
  colnames(X) <- c(outer(1:dim(RUN)[1],1:dim(RUN)[2],\(i,j) paste0("(",i,",",j,")")))
  return(X)
}

X <- rbind(collapse_grid(run1),
           collapse_grid(run2),
           collapse_grid(run3),
           collapse_grid(run4))
# rm(list = c("run1","run2","run3","run4"))

par(mfrow=c(1,3), mar=c(5,2,4,2)+0.1, oma=c(0,4,0,0))  
## minimum of all cells (Final #1)
extreme.level.minimum <- tapply(apply(X,1,min), rep(month.index, times = 165*4), 
                                \(vec) quantile(vec, c(0.9,0.99,0.999,1)))
extreme.level.minimum <- do.call("rbind", extreme.level.minimum)
plot(1:12, extreme.level.minimum[,1], type = "l", ylim = c(0,1.7), lwd = 2, 
     main = "(a) Minimum of all 25 cells", cex.main = 1.5,
     xaxt = "n", xlab = "Month", ylab = "", cex.lab=1.5, col="red")
axis(1, 1:12, 1:12, las=1)
matlines(1:12, extreme.level.minimum[,-1], col=c("blue","green","magenta"), lty = 1, lwd = 2)
abline(h = 1.7); text(3,1.7,labels = "1.7 Leadbetters", pos = 1, cex = 1.5)
abline(v = c(4.5,10.5), lty=2)
# legend("topright", legend = c("max", expression(tau == 0.999), expression(tau == 0.99), expression(tau == 0.9)), 
#        lwd = 2, col = c("red", "blue", "green","magenta")[4:1], cex=1.2)

## X_(6) of all cells (Final #2)
extreme.level.sixth <- tapply(apply(X,1,\(vec) sort(vec)[20]), rep(month.index, times = 165*4), 
                              \(vec) quantile(vec, c(0.9,0.99,0.999,1)))
extreme.level.sixth <- do.call("rbind", extreme.level.sixth)
plot(1:12, extreme.level.sixth[,1], type = "l", ylim = c(0,5.7), lwd = 2, 
     main = "(b) Sixth largest among 25 cells", cex.main = 1.5,
     xaxt = "n", xlab = "Month", ylab = "", cex.lab=1.5, col="red")
axis(1, 1:12, 1:12, las=1)
matlines(1:12, extreme.level.sixth[,-1], col=c("blue","green","magenta"), lty = 1, lwd = 2)
abline(h = 5.7); text(3,5.7,labels = "5.7 Leadbetters", pos = 1, cex = 1.5)
abline(v = c(4.5,10.5), lty=2)
# legend("topright", legend = c("max", expression(tau == 0.999), expression(tau == 0.99), expression(tau == 0.9)), 
#        lwd = 2, col = c("red", "blue", "green","magenta")[4:1], cex=1.2)

## X_(3) of all cells (Final #3)
extreme.level.third <- tapply(apply(X,1,\(vec) sort(vec)[23]), rep(month.index, times = 165*4), 
                              \(vec) quantile(vec, c(0.9,0.99,0.999,1)))
extreme.level.third <- do.call("rbind", extreme.level.third)
plot(1:12, extreme.level.third[,1], type = "l", ylim = c(0,8), lwd = 2, 
     main = "(c) Third largest among 25 cells", cex.main = 1.5,
     xaxt = "n", xlab = "Month", ylab = "", cex.lab=1.5, col="red")
axis(1, 1:12, 1:12, las=1)
matlines(1:12, extreme.level.third[,-1], col=c("blue","green","magenta"), lty = 1, lwd = 2)
abline(h = 5); text(3,5,labels = "5 Leadbetters", pos = 1, cex = 1.5)
abline(v = c(4.5,10.5), lty=2)
legend("topright", legend = c("max", expression(tau == 0.999), expression(tau == 0.99), expression(tau == 0.9)), 
       lwd = 2, col = c("red", "blue", "green","magenta")[4:1], cex=1.5)
mtext("Monthly quantiles (Leadbetters)", side=2, outer=TRUE, line=1)




## temporal tail dependence check ## 

# lag = 1
res_chi <- chiplot(cbind(run1[1,1,-1], run1[1,1,-60225]))

df_chi <- tibble(
  quantile = res_chi$quantile,
  chi = res_chi$chi[, 2],
  chilow = res_chi$chi[, 1],
  chiupp = res_chi$chi[, 3]
)

# Plot for chi
p_chi <- ggplot(df_chi, aes(x = df_chi$quantile, y = df_chi$chi)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = df_chi$chilow, ymax = df_chi$chiupp), alpha = 0.2, fill = "gray30") +
  xlim(0.5, 1) +
  labs(title = "(a) Time lag of one day", x = expression(u), y = expression(chi[1](u))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), 
    plot.title = element_text(face = "bold", size = 20)
  )


# lag = 2
res_chi2 <- chiplot(cbind(run1[1,1,-c(1:2)], run1[1,1,-c(60224:60225)]))


df_chi2 <- tibble(
  quantile = res_chi2$quantile,
  chi = res_chi2$chi[, 2],
  chilow = res_chi2$chi[, 1],
  chiupp = res_chi2$chi[, 3]
)

# Plot for chi
p_chi2 <- ggplot(df_chi2, aes(x = df_chi2$quantile, y = df_chi2$chi)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = df_chi2$chilow, ymax = df_chi2$chiupp), alpha = 0.2, fill = "gray30") +
  xlim(0.5, 1) +
  labs(title = "(b) Time lag of two days", x = expression(u), y = expression(chi[2](u))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 20)
  )


grid.arrange(p_chi, p_chi2, ncol=2)


## spatial tail dependence check ## 

# (1,1) vs (1,2)
res_chi_sp <- chiplot(cbind(run1[1,1,], run1[1,2,]))

df_chi_sp <- tibble(
  quantile = res_chi_sp$quantile,
  chi = res_chi_sp$chi[, 2],
  chilow = res_chi_sp$chi[, 1],
  chiupp = res_chi_sp$chi[, 3]
)

# Plot for chi
p_chi_sp <- ggplot(df_chi_sp, aes(x = df_chi_sp$quantile, y = df_chi_sp$chi)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = df_chi_sp$chilow, ymax = df_chi_sp$chiupp), alpha = 0.2, fill = "gray30") +
  xlim(0.5, 1) +
  labs(title = "(b) Cells (1,1) and (1,2)", x = expression(u), y = expression(chi(u))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), 
    plot.title = element_text(face = "bold", size = 15, hjust=0.5)
  )


# (1,1) vs (3,3)
res_chi_sp2 <- chiplot(cbind(run1[1,1,], run1[3,3,]))

df_chi_sp2 <- tibble(
  quantile = res_chi_sp2$quantile,
  chi = res_chi_sp2$chi[, 2],
  chilow = res_chi_sp2$chi[, 1],
  chiupp = res_chi_sp2$chi[, 3]
)

# Plot for chi
p_chi_sp2 <- ggplot(df_chi_sp2, aes(x = df_chi_sp2$quantile, y = df_chi_sp2$chi)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = df_chi_sp2$chilow, ymax = df_chi_sp2$chiupp), alpha = 0.2, fill = "gray30") +
  xlim(0.5, 1) +
  labs(title = "(c) Cells (1,1) and (3,3)", x = expression(u), y = expression(chi(u))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), 
    plot.title = element_text(face = "bold", size = 15, hjust=0.5)
  )



grid.arrange(p_chi, p_chi_sp, p_chi_sp2, ncol=3)


