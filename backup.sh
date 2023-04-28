#!/usr/bin/env sh

git push -u origin master
rclone copy -v . box_ucdavis:ClusterSearch


