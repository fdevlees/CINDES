
import sqlite3 as lite

dbname = 'omega_stab.db'

class DB(object):
    def __init__(self, dbname):
        try:
            self.con = lite.connect(dbname)
        except lite.Error as e:
            print("error opening db")
            raise
        else:
            self.cur = self.con.cursor()
        return

    def __enter__(self):
        print("enter db...")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.con:
            print("exiting db...")
            self.con.commit()
            self.cur.close()
            self.con.close()

    def get_tablenames(self):
        query = 'select name from sqlite_master where type="table";'
        names, _ = self.do_query(query)
        return [ item[0] for item in names ]

    def get_table(self,tablename):
        query = 'select * from {}'.format(tablename)
        table, description = self.do_query(query)
        return [description] + [ item for item in table ]

    def do_query(self,sql):
        result = self.cur.execute(sql)
        description = list([x[0] for x in result.description])
        return result.fetchall(), description

def to_tablebin( table ):
    import pickle
    ftable = []
    for item in table[1:]:
        index = '_'.join( item[2:7] )
        value = float(item[7])
        ftable.append( [ index, value ] )
    print("len table:", len(ftable))

    with open('tablebin','wb') as f:
        pickle.dump( ftable, f)

def main():
    with DB(dbname) as mydb:
        print("tables:", mydb.get_tablenames())
        names = mydb.get_tablenames()
        for i,name in enumerate(names):
            print("i:", i, "name:", name)
            table1 = mydb.get_table( name )
            for item in table1[:3]:
                print(" | ".join(map(str,item)))
            print("n:", len(table1))
            print()

        to_tablebin( mydb.get_table( names[2] ) )





if __name__=="__main__":
    main()

