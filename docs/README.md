# Data Preprocessing

The flowchart below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The coloured silos, rounded boxes, and grey boxes, represent data, R scripts, and directories respectively. The data is coloured according to the [medallion architecture](https://www.databricks.com/glossary/medallion-architecture), where bronze represents raw data, silver represents cleansed and processed data, and gold represents data which are consumption ready (i.e. directly used within statistical models).

```mermaid
flowchart LR 
	subgraph data
		d11[(FlowR15_R18_VoIs_221118.csv)]
		d12[(flowR19_VOIs.dta)]
 
        d21[(Quest_R015_R020_VOIs_May062022.csv)]
        d22[(Allpcr_data_for_R015_R020_study_ids.xlsx)]
        d23[(R016_R020_Data_for_HIVnegatives.csv)]
        d24[(RCCS_1strndquest_participation.dta)]


	end 
	style d11 fill:#CD7F32
	style d12 fill:#CD7F32
	style d21 fill:#CD7F32
	style d22 fill:#CD7F32
	style d23 fill:#CD7F32
	style d24 fill:#CD7F32
	
	subgraph src
		s1(get_census_eligible_count.R)
        s2(process_data.R)
        s3(get_participation_rates.R) 
	end
	
	subgraph outputs
		o11[(census_eligible_individuals_230514.csv)]
        o21[(all_participants_hivstatus_vl_230502.csv)]
        o31[(participation_rates_230517.rds)]

		style o11 fill:#F5F5F5
		style o21 fill:#F5F5F5
		style o31 fill:#F5F5F5
		style o11 color:#000000
		style o21 color:#000000
		style o31 color:#000000
	end
	
	d11 -- load --> s1
	d12 -- load --> s1
	s1 -- create --> o11
  
	d21 -- load --> s2
	d22 -- load --> s2
	d23 -- load --> s2
	d24 -- load --> s2
	s2 -- create --> o21
 
	o11 -- load --> s3
	o21 -- load --> s3
	s3 -- create --> o31



```


Successively, the outputs of the preprocessing are used in order to run our models.
(To finish)

```mermaid
flowchart LR 
	subgraph data
		o11[(census_eligible_individuals_230514.csv)]
        o21[(all_participants_hivstatus_vl_230502.csv)]

    style o11 fill:#F5F5F5
    style o21 fill:#F5F5F5
    style o11 color:#000000
    style o21 color:#000000
	end 
 
	subgraph scripts
		sc1(VL_run_stan.R)
        sc2(VL_postprocessing.R)
	end
	
	subgraph results
		r11[(TODO)]
        r21[(TODO)]

		style r11 fill:#ffd700
		style r21 fill:#ffd700
		style r11 color:#ffffff
		style r21 color:#ffffff
	end
	
	o11 -- load --> sc1
	sc1 -- create --> r11
  
	o21 -- load --> sc2
	sc2 -- create --> r21
```
