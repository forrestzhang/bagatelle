from pymongo import MongoClient
import json

client = MongoClient('localhost', 27017)
db = client.pubmed
collection = db.rice

with open("/Users/forrest/Downloads/pubmed_result.json") as inio:
    for i in inio:
        i = i.rstrip()
        ij = json.loads(i)

        try:
            result = collection.insert_many(ij)
            result.inserted_ids

        except:

            pass
