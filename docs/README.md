# Data Preprocessing

The flowchart below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The coloured silos, rounded boxes, and grey boxes, represent data, R scripts, and directories respectively. The data is coloured according to the [medallion architecture](https://www.databricks.com/glossary/medallion-architecture), where bronze represents raw data, silver represents cleansed and processed data, and gold represents data which are consumption ready (i.e. directly used within statistical models).

```mermaid
flowchart LR 
	subgraph data
		d1[(FlowR15_R18_VoIs_221118.csv)]
		d2[(flowR19_VOIs.dta)]
	end 
	style d1 fill:#CD7F32
	style d2 fill:#CD7F32
	
	subgraph src
		s1(get_census_eligible_count.R)
	end
	
	subgraph outputs
		o1[(census_eligible_individuals_230514.csv)]
		style o1 fill:#FFD700
		style o1 fill:#F5F5F5
	end
	
	d1 -- load --> s1
	d2 -- load --> s1
	s1 -- create --> o1
```



