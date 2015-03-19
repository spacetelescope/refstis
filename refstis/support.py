import sys, os  #fix this sometime  Remove this and run Dark_Monitor.py and see what is wrong

#---------------------------------------------------------------------------
#---------------------SQL
#---------------------------------------------------------------------------

def createXmlFile(ftp_dir=None,set=None, file_type="all", archive_user=None, archive_pwd=None, email=None, host='science3.stsci.edu', ftp_user=None, ftp_pwd=None):

    import sys, traceback, os, string, errno, glob, httplib, urllib, time

    exposure_template = string.Template('\
<?xml version=\"1.0\"?> \n \
<!DOCTYPE distributionRequest SYSTEM \"http://archive.stsci.edu/ops/distribution.dtd\"> \n \
  <distributionRequest> \n \
    <head> \n \
      <requester userId = \"$archive_user\" email = \"$email\" archivePassword = \"$archive_pwd" source = "starview"/> \n \
      <delivery> \n \
        <ftp hostName = \"$host\" loginName = \"$ftp_user\" loginPassword = \"$ftp_pwd\" directory = \"$ftp_dir\" secure = \"true\" /> \n \
      </delivery>    <process compression = \"none\"/> \n \
    </head> \n \
    <body> \n \
      <include> \n \
        <select> \n \
          $types \n \
        </select> \n \
        $datasets \n \
      </include> \n \
    </body> \n \
  </distributionRequest> \n' )

    print archive_user,archive_pwd,ftp_user,ftp_pwd,ftp_dir,email

    if not archive_user: archive_user = raw_input('I need your archive username: ')
    if not archive_pwd: archive_pwd = raw_input('I need your archive password: ')
    if not ftp_user: ftp_user = raw_input('I need the ftp username: ')
    if not ftp_pwd: ftp_pwd = raw_input('I need the ftp password: ')
    if not ftp_dir: ftp_dir = raw_input('Where would you like the files delivered: ')
    if not email: email = raw_input('Please enter the email address for the notification: ')
    datasets = "\n"

    if file_type == "all":
        types = "<suffix name=\"*\" />"
    elif file_type == "caldq":
        types = "<retrievalKeyword name=\"DataQuality\" /> \n           <retrievalKeyword name=\"JitterFiles\" />"
    elif file_type == "lpf":
        types = "<suffix name=\"LPF\" /> \n           <suffix name=\"SPT\" />"
    elif file_type == "raw":
        types = "<retrievalKeyword name=\"Uncalibrated\" />"
    elif file_type == "raw+dq+jit":
        types = "<retrievalKeyword name=\"Uncalibrated\" />\n           <retrievalKeyword name=\"JitterFiles\" /> \n           <retrievalKeyword name=\"DataQuality\" /> \n           <suffix name=\"MMD\" />\n           <suffix name=\"SPT\" />"
    elif file_type == "crj":
        types = "<suffix name=\"CRJ\" />\n"
    else:
        types = "<suffix name=\"*\" />\n"

    request_template = exposure_template
    request_str = request_template.safe_substitute( archive_user = archive_user, \
                                                        archive_pwd = archive_pwd, \
                                                        email = email, \
                                                        host = host, \
                                                        ftp_dir = ftp_dir, \
                                                        ftp_user = ftp_user, \
                                                        ftp_pwd = ftp_pwd, \
                                                        types = types)
    request_template = string.Template(request_str)
    for item in set:
        datasets = datasets+"          <rootname>%s</rootname>\n" % (item)
    xml_file = request_template.safe_substitute(datasets = datasets)
    return xml_file

def submitXmlFile(xml_file,server):
    import sys, traceback, os, string, errno, httplib, urllib

    if server == "dmsops1.stsci.edu":
        params = urllib.urlencode({'dadshost':'dmsops1.stsci.edu','dadsport': 4703, 'request': xml_file})
    else:
        params = urllib.urlencode({'dadshost':'sthubbins.stsci.edu','dadsport': 4703, 'request': xml_file})

    headers = {"Accept": "text/html", "User-Agent":"%sDADSAll" % os.environ.get("USER")}
    req = httplib.HTTPSConnection("archive.stsci.edu")
    req.request("POST", "/cgi-bin/dads.cgi", params, headers)
    resp = req.getresponse()
    textlines = resp.read()
    print 'Request Submitted'
    return textlines

class SybaseInterface:
        import sys
	def __init__(self,server,database,user="york",fields="",tables="",where="",order="",distinct=False):
                import sys
		self.server = server
		self.dbname = database
		self.fields_txt = fields
		self.tables_txt = tables
		self.where_txt = where
		self.order_txt = order
		self.query_result = []
		self.query_dict = {}
		self.distinct = distinct
		self.user = user

	#find these on a solaris machine in /usr/local/sybase/stbin/{username}.dat
	passwords = {"ROBBIE": "janawtmaxi", "ZEPPO": "omwaddgazeba" , "HARPO": "W@@s6L9Np@ub"}

	def findLongest(self,name,strings):
		theLongest = len(name)
		for string in strings:
			if len(string) > theLongest:
				theLongest = len(string)
		return theLongest

	def csvPrint(self,result=None):
		if result == None:
			result = self.resultAsDict()
		columns = result["COLUMNS"]
		str = columns[0]
		for name in columns:
			str += "," + name
		str += "\n"
		for i,name in enumerate(result[columns[0]]):
			str += result[columns[0]][i]
			for cname in columns[1:]:
				str += "," + result[cname][i]
			str += "\n"
		return str

	def prettyPrint(self,result=None):
		if result == None:
			result = self.resultAsDict()
		columns = result["COLUMNS"]
		longest_values = []
		str = "|"
		for name in columns:
			longest_values.append(self.findLongest(name,result[name]))
			str = str + "%-*s|" % (longest_values[-1],name)
		str = str + "\n"
		for i,name in enumerate(result[columns[0]]):
			str = str + "|"
			for j,cname in enumerate(columns):
				str = str + "%-*s|" % (longest_values[j],result[cname][i])
			str = str + "\n"
		return str

	def setDistinct(self,distinct):
		self.distinct = distinct

	def setFieldsText(self,fields):
		self.fields_txt = fields

	def setFields(self,fields):
		self.fields_txt = fields[0]
		for field in fields[1:]:
			self.fields_txt = "%s,%s" % (self.fields_txt,field)

	def setTablesText(self,tables):
		self.tables_txt = tables

	def setTables(self,tables):
		self.tables_txt = "%s" % (self.convertTables(tables))

	def setOrder(self,order):
		self.order_txt = "%s" % (order[0])

	def convertTables(self,item):
		if item[0] == "TABLE":
			text_str = item[1]
		else:
			operation = item[0]
			table_a = self.convertTables(item[1])
			table_b = self.convertTables(item[2])
			on_list = self.getList(item[3])
			text_str = "%s %s %s ON %s" % (table_a,operation,table_b,on_list)
		return text_str

	def getList(self,item):
		text_str = "%s=%s" % (item[0][0],item[0][1])
		for pair in item[1:]:
			text_str = text_str + " AND %s=%s" % (pair[0],pair[1])
		return text_str

	def setWhereText(self,where):
		self.where_txt = where

	def setWhere(self,where):
		self.where_txt = "(%s)" % (self.convertWhere(where))

	def convertWhere(self,item):
		if item[0] == "ATOM":
			keyword = item[1]
			test = item[2]
			if test != "BETWEEN":
				if item[3] == 'T':
					value = "'%s'" % (item[4])
				elif item[3] == 'D':
					value = "%d" % (item[4])
				elif item[3] == 'F':
					value = "%e" % (item[4])
				else:
					value = "%s" % (item[4])
			else:
				if item[3] == 'D':
					value = "%d AND %d" % (item[4],item[5])
				elif item[3] == 'F':
					value = "%e AND %e" % (item[4],item[5])
			text_str = "%s %s %s" % (keyword,test,value)
		else:
			join_str = item[0]
			text_str = "(%s)" % (self.convertWhere(item[-1]))
			for i in range(len(item)-2):
				j = len(item) - i - 2
				text_str = "(%s) %s %s" % (self.convertWhere(item[j]),join_str,text_str)
		return text_str

	def resultAsText(self):
		if sys.platform == "darwin" or sys.platform == "linux2":
			result_txt = self.prettyPrint()
		elif sys.platform == "sunos5":
			result_txt = ""
			for line in self.query_result:
				result_txt = "%s%s\n" % (result_txt,line)
		else:
			result_txt = self.prettyPrint()
		return result_txt

	def resultAsDict(self):
		names = []
                if sys.platform == "darwin" or sys.platform == "linux2":
			for item in self.query_result[2].split("|||"):
				names.append(item.strip().replace('1> 2> ',''))
                elif sys.platform == "sunos5":
			for item in self.query_result[0].split("|")[1:-1]:
				names.append(item.strip())
			lengths = []
			for item in self.query_result[1].split("|")[1:-1]:
				lengths.append(len(item))
		mappings = {}#{"COLUMNS":names}
		for name in names:
			mappings[name] = []
		if sys.platform == "darwin" or sys.platform == "linux2":
			for line in self.query_result[3:]:
				items = line.split("|||")
				for name,item in zip(names,items):
                                        mappings[name].append(item)
		elif sys.platform == "sunos5":
			for i in range(len(self.query_result[2:])):
				j = i + 2
				items = self.lineSplit(self.query_result[j],lengths)
                                for k in range(len(names)):
					mappings[names[k]].append(items[k])
		return mappings

	def lineSplit(self,line,lengths):
		counter = 1
		items = []
		for length in lengths:
			items.append(line[counter:(counter+length)])
			counter = counter + length + 1
		for i in range(len(items)):
			items[i] = items[i].strip()
		return items

	def buildString(self):
		query_string = "SELECT "
		if self.distinct:
			query_string = query_string + "DISTINCT "
		query_string = query_string + "%s FROM %s" % (self.fields_txt,self.tables_txt)
		if self.where_txt != "":
			query_string = query_string + " WHERE %s" % (self.where_txt)
		if self.order_txt != "":
			query_string = query_string + " ORDER BY %s" % (self.order_txt)
		return query_string

	def doQuery(self,query=""):
                if not sys.platform.startswith('linux'):
                        errmsg='New data can only be retrieved from a linux system.\n'
                        errmsg+='\n'
                        errmsg+='Please move to a linux systen or get the data manually and run again \n'
                        errmsg+='with the --no_collect option.'
                        sys.exit(errmsg)

		if query:
			query_string = query
		else:
			query_string = self.buildString()
#		print query_string
		if sys.platform == "darwin" or sys.platform == "linux2":
			query_string = query_string + "\ngo\n"
			transmit,receive,error = os.popen3("tsql -S%s -D%s -U%s -P%s -t'|||'" % (self.server,self.dbname,self.user,self.passwords[self.server]))
		elif sys.platform == "sunos5":
			query_string = query_string + "\ngo\n"
			transmit,receive,error = os.popen3("/usr/local/syb_12.5/OCS-12_5/bin/isql -S%s -D%s -w50000 -s'|'" % (self.server,self.dbname))
		else:
			query_string = query_string + ";\n"
			transmit,receive,error = os.popen3("isql -S%s -D%s -s'|||' -w50000 -b" % (self.server,self.dbname))
		transmit.write(query_string)
		transmit.close()
		self.query_result = receive.readlines()
		receive.close()
		self.error_report = error.readlines()
		error.close()
#		print self.query_result
		for i in range(len(self.query_result)):
			self.query_result[i] = self.query_result[i].strip()
		if sys.platform == "darwin" or sys.platform == "linux2":
			self.query_result = self.query_result[4:-1]
			if sys.platform == "darwin":
				self.query_result = self.query_result[2:]
#			print self.query_result
			if "affected" in self.query_result[-1]:
				self.query_result = self.query_result[:-1]
			self.query_result[0] = self.query_result[0][6:]
		elif sys.platform == "sunos5":
			self.query_result = self.query_result[:-2]

#---------------------------------------------------------------------------

def decimal_year(month,day,year):
    '''
    Returns decimal day for given mm,dd,yyyy
    '''
    day=int(day)
    month=int(month)
    year=int(year)
    if month==1: decimal_day=day+0
    elif month==2: decimal_day=day+31
    elif month==3: decimal_day=day+59
    elif month==4: decimal_day=day+90
    elif month==5: decimal_day=day+120
    elif month==6: decimal_day=day+151
    elif month==7: decimal_day=day+181
    elif month==8: decimal_day=day+212
    elif month==9: decimal_day=day+243
    elif month==10: decimal_day=day+273
    elif month==11: decimal_day=day+304
    elif month==12: decimal_day=day+334
    else:
        print 'Input Error'
    decimal_day=decimal_day/365.0
    decimal_date=year+decimal_day
    return decimal_date

#---------------------------------------------------------------------------

def sigma_clip(input_array, sigma=5, iterations=30, print_message=False):
    import os, numpy
    try: median
    except NameError: from numpy import where,median

    input_array=numpy.array(input_array)
    if len(input_array.shape)>1:
    	input_array=input_array.flatten()
    for iteration in xrange(iterations):
    	if print_message:
    	    os.write(1,'\b\b\b\b\b\b\b\b\bPass: %d'%(iteration))

	std=input_array.std()
	ceiling=numpy.median(input_array)+sigma*std
	floor=numpy.median(input_array)-sigma*std
        index=numpy.where((input_array<ceiling) & (input_array>floor))[0]

        if len(index)==len(input_array):
            #print 'Done'
            break
        else:
	    input_array=input_array[index]

    if print_message:
        print '\n'
    return numpy.median(input_array),input_array.mean(),input_array.std()

#---------------------------------------------------------------------------

def mjd_to_greg(mjd):
   #This comes from http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
   JD = mjd + 2400000.5
   JD = int(JD)
   L= JD+68569
   N= 4*L/146097
   L= L-(146097*N+3)/4
   I= 4000*(L+1)/1461001
   L= L-1461*I/4+31
   J= 80*L/2447
   K= L-2447*J/80
   L= J/11
   J= J+2-12*L
   I= 100*(N-49)+I+L
   Year = I
   Month = J
   Day = K
   month_to_day = {'0': 0,'1':31, '2':59, '3':90, '4':120, '5':151, '6':181, '7':212, '8':243, '9':273, '10':304, '11':334, '12':365}
   tot_day = (month_to_day[str(int(Month)-1)] + Day)
   day_in_year = 365.0
   if (Month >= 2) & (Year%4 == 0): #for leap year
      tot_day = tot_day + 1.0
      day_in_year = 365.0 + 1.0
   frac_day = tot_day / day_in_year
   fractional_year = Year + frac_day
   return (Year, Month, Day, fractional_year)

#---------------------------------------------------------------------------

class Logger(object):
    def __init__(self,filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()

#---------------------------------------------------------------------------

def send_email(subject=None,message=None,from_addr=None,to_addr=None):
    '''
    Send am email via SMTP server.
    This will not prompt for login if you are alread on the internal network.
    '''
    import os
    import getpass
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    users_email=getpass.getuser()+'@stsci.edu'

    if not subject:
	subject='Message from %s'%(__file__)
    if not message:
	message='You forgot to put a message into me'
    if not from_addr:
        from_addr=users_email
    if not to_addr:
        to_addr=users_email

    svr_addr='smtp.stsci.edu'
    msg = MIMEMultipart()
    msg['Subject']=subject
    msg['From']=from_addr
    msg['To']=to_addr
    msg.attach(MIMEText(message))
    s = smtplib.SMTP(svr_addr)
    s.sendmail(from_addr, to_addr, msg.as_string())
    s.quit()
    print '\nEmail sent to %s \n' %(from_addr)
