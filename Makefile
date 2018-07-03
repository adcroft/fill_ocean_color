FILES = S19980012010031.L3m_MC_CHL_chlor_a_9km.nc \
        S19980322010059.L3m_MC_CHL_chlor_a_9km.nc \
        S19980602010090.L3m_MC_CHL_chlor_a_9km.nc \
        S19980912010120.L3m_MC_CHL_chlor_a_9km.nc \
        S19981212010151.L3m_MC_CHL_chlor_a_9km.nc \
        S19981522010181.L3m_MC_CHL_chlor_a_9km.nc \
        S19981822010212.L3m_MC_CHL_chlor_a_9km.nc \
        S19982132010243.L3m_MC_CHL_chlor_a_9km.nc \
        S19972442010273.L3m_MC_CHL_chlor_a_9km.nc \
        S19972742010304.L3m_MC_CHL_chlor_a_9km.nc \
        S19973052010334.L3m_MC_CHL_chlor_a_9km.nc \
        S19973352010365.L3m_MC_CHL_chlor_a_9km.nc
TOPO = GEBCO_08_v1.nc
TARG = seawifs-clim-1997-2010.nc
HASH = hash.md5

# Main target checks the checksums
check: $(HASH)
	md5sum -c $(HASH)
Check: check
	md5sum -c raw.md5
	md5sum -c gebco.md5

# Fetch data from web
S%.nc:
	wget https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/$@

# Fetch topography file from archive (not easily available on web)
$(TOPO):
	cp /archive/gold/datasets/topography/$@ .
	md5sum -c gebco.md5

# Fill and join data
$(TARG): $(TOPO) $(FILES)
	./fill_and_join_chlor_a.py $@ $(TOPO) $(FILES)

# Record checksums
$(HASH): | $(TARG)
	md5sum $(TARG) > $@
raw.md5: $(FILES)
	md5sum $(FILES) > $@

# Clean up
clean:
	rm *.nc
