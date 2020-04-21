#!/bin/bash
function SetUpServer(){
git fetch origin master
git reset --hard origin/master
}
SetUpServer 
