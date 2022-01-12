out <- prospect5(1.4, 40, 10, 0, 0.01, 0.01)
plot(c(400, 2500), c(0, 1), type = "n",
     xlab = "Wavelength (nm)",
     ylab = "Spectra")
lines(400:2500, out$reflectance, col = "red")
lines(400:2500, 1 - out$transmittance, col = "blue")
legend("top", legend = c("reflectance", "1 - transmittance"),
       lty = 1, col = c("red", "blue"))