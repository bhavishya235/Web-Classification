#!/bin/bash

ff=0
count=0
`mkdir AllVector2/$ff`
for i in `ls -1 AllVector`
	do
		let count=count+1
		if [ $count -eq 100 ]
		then
			let ff=ff+1
			mkdir AllVector2/$ff
			let count=0
		fi
		cp AllVector/$i AllVector2/$ff/$i
	done
