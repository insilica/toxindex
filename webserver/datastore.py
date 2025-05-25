import psycopg2
import psycopg2.extras
import os
import logging

# Read environment variables for PostgreSQL connection
dbhost = os.getenv('DB_HOST')
dbport = os.getenv('DB_PORT')
dbname = os.getenv('DB_NAME')
dbuser = os.getenv('DB_USER')
dbpass = os.getenv('DB_PASSWORD')

def get_connection():
    con = psycopg2.connect(host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)
    psycopg2.extras.register_uuid(conn_or_curs=con)
    return con

def find(query, param=None):
    res = None
    con = None
    cur = None
    try:
        con = get_connection()
        cur = con.cursor(cursor_factory=psycopg2.extras.DictCursor) # Use DictCursor here
        if param is None:
            cur.execute(query)
        else:
            cur.execute(query, param)
        res = cur.fetchone()
    except Exception as e:
        logging.debug("error executing query")
        logging.debug(query)
        logging.debug(e)
    finally:
        if cur: cur.close()
        if con: con.close()
        return res

def execute(query, param=None):
    con = None
    try:
        con = get_connection()
        cur = con.cursor()
        if param is None:
            cur.execute(query)
        else:
            cur.execute(query, param)
        con.commit()
    except Exception as e:
        logging.debug("error executing query")
        logging.debug(query)
        logging.debug(e)
    finally:
        if con: con.close()
        logging.debug("done")

def find_all(query, param=None):
    res = []
    con = None
    cur = None
    try:
        con = get_connection()
        cur = con.cursor(cursor_factory=psycopg2.extras.DictCursor)  # Use DictCursor here
        if param is None:
            cur.execute(query)
        else:
            cur.execute(query, param)
        res = cur.fetchall()
    except Exception as e:
        logging.debug("error executing query")
        logging.debug(query)
        logging.debug(e)
    finally:
        if cur: cur.close()
        if con: con.close()
        return res