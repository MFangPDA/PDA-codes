%% to load NC data
function data = subf_readNC(nc_path,var_name)
tas_ncid=netcdf.open(nc_path,'NC_NOWRITE');
tas_varid=netcdf.inqVarID(tas_ncid,var_name);  
data=netcdf.getVar(tas_ncid,tas_varid);        
netcdf.close(tas_ncid);
end