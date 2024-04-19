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
from django.urls import path, re_path
from drf_yasg.views import get_schema_view
from drf_yasg import openapi
from rest_framework import permissions


from blast import views
from blast.open_api.MyView import MyView
from blast.open_api.alignement.AlignementView import AlignementView

schema_view = get_schema_view(
   openapi.Info(
      title="Your Project API",
      default_version='v1',
      description="API documentation for Your Project",
   ),
   public=True,
   permission_classes=(permissions.AllowAny,),
)
urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.blast),
    path('blast', views.blast),
    path('documentation', views.documentation),
    re_path(r'^swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    path('api_decumentation/', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    path('redoc/', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),
    path('alignement_viewer', views.alignement_viewer),
    path('submit_sequence', views.submit_sequence),
    path('submit_sequence_query', views.submit_sequence_query),
    path('result_viewer', views.result_viewer),


    path('alignement_viewer', AlignementView.as_view(), name='alignement_viewer'),


    path('hello_world/', MyView.as_view(), name='hello_world'),
]
