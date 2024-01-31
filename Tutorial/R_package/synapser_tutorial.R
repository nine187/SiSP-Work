library(synapser) 
library(synapserutils) 
## 
## TERMS OF USE NOTICE:
##   When using Synapse, remember that the terms and conditions of use require that you:
##   1) Attribute data contributors when discussing these data or results from these data.
##   2) Not discriminate, identify, or recontact individuals or groups represented by the data.
##   3) Use and contribute only data de-identified to HIPAA standards.
##   4) Redistribute data only under these same terms of use.
synLogin(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcwMzEyODkyNiwiaWF0IjoxNzAzMTI4OTI2LCJqdGkiOiI0NjM5Iiwic3ViIjoiMzQ4OTc2MCJ9.gawO8XMu4riidhSTqSqGlW17mOcdz4SNc8hH-jhp87n1DB29ABehKjTW-u5RE3JjZt0nhLrSMMk2bcJ5Z6mG--7Nzvxwf0lhtYqS1zb2GtKct2_k8hf2uuD9lGi1krL8ugf06rVGy_M-kY5PiIH09LC98llISdTWPbyvIiDGQJtAfgbEJHk7_7bzU1HNa5rO1zwaNf4ZXxsLf-S0NpdxFYfp3xbrc1V8WsvV9YWoqpmfjISXZY5nJNwsQsWByx7cfPJuuh1XFEy4i-KNhUWe3H6jnkvFjy8PoMkKGE3XO7c-O2aAicC4lnYZ2mqLy4kq--qHaqSB0fQjEWvblPmrsg")
## NULL
#get the authenthication token from https://www.synapse.org/#!PersonalAccessTokens:

files <- synapserutils::syncFromSynapse('syn2634724') 