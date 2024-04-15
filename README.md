# bioBlastAutomation :rocket:

`bioBlastAutomation` is a Python-based project designed to automate the process of performing BLAST (Basic Local Alignment Search Tool) operations. BLAST is a bioinformatics algorithm used to compare primary biological sequence information, such as the amino-acid sequences of proteins or the nucleotides of DNA and/or RNA sequences.

This project utilizes the Biopython library, a set of tools for biological computation, and django, a web framework for Python. It also employs Plotly for data visualization.

The `bioBlastAutomation` project includes utilities for validating and identifying the type of sequences (DNA, RNA, Protein) and for checking the compatibility of the sequence with the BLAST program. It also includes functionality for saving BLAST results in different formats (XML, JSON, CSV).

A unique feature of this project is its web interface, which displays a bar chart of sequence alignments. These alignments are read from a file and visualized using Plotly.

The project is managed with pip and has dependencies on several libraries, as listed in the `requirements.txt` file.



## Installation :wrench:

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

```bash
# Clone the repository
git clone https://github.com/abdelhak-zaaim/BioBlastAutomation.git

# Navigate to the project directory
cd BioBlastAutomation

# Install the required packages
pip install -r requirements.txt

```

## RabbitMQ Installation and Configuration

RabbitMQ is an open-source message broker that we use in this project for task queueing with Celery. Here's how you can install and configure it:

### Installing RabbitMQ

1. Install Homebrew if you haven't already. You can install it by running the following command in your terminal:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

```

2. Install RabbitMQ by running the following command in your terminal:

```bash
brew install rabbitmq
```
3. add the following line to your `.bashrc` or `.zshrc` file:

```bash
export PATH=$PATH:/usr/local/sbin
source ~/.zshrc
```

4. open the Constants.py file and change the value of the `RABBITMQ_USERNAME` and `RABBITMQ_PASSWORD` variables to your RabbitMQ username and password.

```python
class Constants:
    MY_CELERY_APP = 'CELERY_APP' # Celery app name could be any name
    MY_CELERY_BROKER = 'amqp://guest:guest@localhost:5672//' # Celery broker URL its default value you may change it

```



5. Start the RabbitMQ server by running the following command in your terminal:

```bash
rabbitmq-server start
```


6. Apply the migrations by executing the commands
    
```bash
python manage.py makemigrations
python manage.py migrate
```
7. Run the application by executing the command (web interface) 
    
```bash
python manage.py runserver
```

## finaly you can access the web interface by visiting the url its shown in the terminal




![Project Image](https://fsdm.zaaim.me/src/images/Screenshot%202024-04-03%20at%2022.02.10.png)

