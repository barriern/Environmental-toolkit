load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

fin = addfile("data/combined_gridded_index_data.nc", "r")
nao = fin->nao
amo = fin->amo
sst = fin->msl
occu = new((/2, dimsizes(nao)/), float)
occu(0, :) = nao
occu(1, :) = amo
occu!0 = "ncl"
occu!1 = "year"

nlat = 80
nlon = 147
ncl = 2
maxlag = 5

printVarSummary(sst)

xleady = esccr(occu(ncl|:, year|:), sst(lat|:, lon|:, year|:), maxlag)
xleady!0="ncl"
xleady!1="lat"
xleady!2="lon"
xleady!3="lag"

yleadx = esccr(sst(lat|:, lon|:, year|:), occu(ncl|:, year|:), maxlag)
yleadx!0="lat"
yleadx!1="lon"
yleadx!2="ncl"
yleadx!3="lag"

output = new((/ncl, nlat, nlon, 2*maxlag+1/), float)
output(:, :, :, 0:maxlag-1) = yleadx(ncl|:, lat|:, lon|:, lag|1:maxlag:-1)
output(:, :, :, maxlag:) = xleady(:, :, :, 0:maxlag)
output&lon = sst&lon
output&lat = sst&lat
lag = ispan(-maxlag, maxlag, 1)
output&lag = lag

fout = addfile("correlation_ncl.nc", "c")
fout->corr = output



end
