from pytube import YouTube
yt = YouTube("https://www.youtube.com/watch?v=PkDr5OogdQo")
yt = yt.get('mp4', '720p')
yt.download('C:/Users/Prasoon Gupta')