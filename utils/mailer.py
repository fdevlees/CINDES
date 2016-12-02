# import smtplib for the actual sending function
import smtplib
# import the email modules needed
from email.mime.text import MIMEText

def sendmail(filename,folder):
    message = 'a job crashed with name: %s in folder %s' %(filename,folder)
    # create a text/plain message
    msg = MIMEText(message)

    me = 'error@hydra.vub.ac.be'
    you = 'jlteunissen@gmail.com'
    msg['Subject'] = 'The contents of %s' %(filename)
    msg['From'] = me
    msg['To'] = you
    # send the message via our own SMTP server
    s = smtplib.SMTP('localhost')
    s.sendmail(me, [you], msg.as_string())
    s.quit()

if __name__ == "__main__":
    filename = 'blabla'
    folder = 'diamant'
    sendmail(filename,folder)

