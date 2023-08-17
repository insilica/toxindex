import psycopg
import os

# Read environment variables for PostgreSQL connection
dbhost = os.environ.get('POSTGRES_HOST')
dbport = os.environ.get('POSTGRES_PORT')
dbname = os.environ.get('POSTGRES_DB')
dbuser = os.environ.get('POSTGRES_USER')
dbpass = os.environ.get('POSTGRES_PASSWORD')

def find(query, param=None):
    res = None
    con = None
    cur = None
    try:
        con = psycopg.connect(host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)
        cur = con.cursor()
        if param is None:
            cur.execute(query)
        else:
            cur.execute(query, param)
        res = cur.fetchone()
    except Exception as e:
        print("error executing query") # TODO actually log this.
        print(query)
        print(e)
    finally:
        if cur: cur.close()
        if con: con.close()
        return res

def execute(query, param=None):
    con = None
    try:
        con = psycopg.connect(host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)
        cur = con.cursor()
        if param is None:
            cur.execute(query)
        else:
            cur.execute(query, param)
        con.commit()
    except Exception as e:
        print("error executing query") # TODO actually log this.
        print(query)
        print(e)
    finally:
        if con: con.close()
        return "done"
