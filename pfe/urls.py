"""
URL configuration for pfe project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.blast, name='blast')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='blast')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path

from blast import views

urlpatterns = [
    path('admin/', admin.site.urls),

    path('blast', views.blast),
    path('documentation', views.documentation),

    path('alignement_viewer', views.alignement_viewer),
    path('submit_sequence', views.submit_sequence),
]

